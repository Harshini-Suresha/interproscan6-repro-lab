// Generates synthetic proteins, runs the miniature scanner, writes data + measured benchmarks.
const fs=require('fs'),crypto=require('crypto');
const PAT=[["PS00018","EF-hand calcium-binding","D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]"],
["PS00028","Zinc finger C2H2","C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H"],
["PS00017","ATP/GTP-binding P-loop","[AG]-x(4)-G-K-[ST]"],
["PS00029","Leucine zipper","L-x(6)-L-x(6)-L-x(6)-L"],
["PS00001","N-glycosylation site","N-{P}-[ST]-{P}"],
["PS00008","N-myristoylation site","G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}"]];
const toRe=s=>new RegExp(s.split("-").map(t=>{const m=t.match(/^x\((\d+)(?:,(\d+))?\)$/);if(m)return"."+"{"+m[1]+(m[2]?","+m[2]:"")+"}";if(t==="x")return".";if(t[0]==="{")return"[^"+t.slice(1,-1)+"]";return t}).join(""),"g");
PAT.forEach(p=>p.push(toRe(p[2])));
function scanOne(id,seq,out){for(const p of PAT){const re=p[3];re.lastIndex=0;let m;while((m=re.exec(seq))){out.push([id,"-",seq.length,"PROSITE-lite",p[0],p[1],m.index+1,m.index+m[0].length,"-","T","-","-","-","-","-"].join("\t"));re.lastIndex=m.index+1}}}
const rng=a=>()=>{a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296};
const AA="ACDEFGHIKLMNPQRSTVWY";
const PL=[["APSTLGKT","PS00017"],["CAACAAALAAAAAAAAHAAAH","PS00028"],["LAAAAAALAAAAAALAAAAAAL","PS00029"],["DKDGDGKIEAAEL","PS00018"],["GKSAGSNGAA",null]];
function synth(n,seed){const r=rng(seed),o=[],truth=[];for(let i=0;i<n;i++){const L=150+Math.floor(r()*450);let a=[];for(let k=0;k<L;k++)a.push(AA[Math.floor(r()*20)]);let s=a.join("");const pl=[];
 for(const [m,sig] of PL){if(r()<.3){const p=Math.floor(r()*(s.length-m.length));s=s.slice(0,p)+m+s.slice(p+m.length);pl.push([p,m,sig])}}
 const id="SYN"+String(i+1).padStart(5,"0");o.push([id,s]);
 for(const [p,m,sig] of pl)if(sig&&s.substr(p,m.length)===m)truth.push([id,sig,p+1,p+m.length].join("\t"))}return{prots:o,truth}}
const fp=l=>crypto.createHash("sha256").update(l.slice().sort().join("\n")).digest("hex").slice(0,12);
function run(P,chunk,order){let list=P.slice();if(order==="reverse")list.reverse();if(order==="shuffle"){const r=rng(99);for(let i=list.length-1;i>0;i--){const j=Math.floor(r()*(i+1));[list[i],list[j]]=[list[j],list[i]]}}
 let peak=process.memoryUsage().rss;const t0=process.hrtime.bigint(),all=[];
 for(let i=0;i<list.length;i+=chunk){for(const p of list.slice(i,i+chunk))scanOne(p[0],p[1],all);peak=Math.max(peak,process.memoryUsage().rss)}
 const ms=Number(process.hrtime.bigint()-t0)/1e6;return{ms,rss:peak/1048576,all}}
const CF=[["chunk1000_in_order",1000,"in"],["chunk100_in_order",100,"in"],["chunk10_shuffled",10,"shuffle"],["chunk500_reversed",500,"reverse"]];
const bench=["dataset,n_proteins,seed,config,repeat,wall_ms,peak_rss_mb,matches,coverage_pct,fingerprint"],truthcmp=["dataset,n_proteins,planted_motifs,recovered,recall_pct,total_matches,chance_matches,coverage_pct"],bysig=["dataset,signature,description,matches,proteins_with_hit"];
const SETS=[[100,1],[500,1],[1000,1],[1000,2],[5000,1],[20000,1]];
for(const [n,seed] of SETS){const name=`synthetic_${n}_seed${seed}`,{prots,truth}=synth(n,seed);
 fs.writeFileSync(`data/${name}.fasta`,prots.map(p=>">"+p[0]+"\n"+p[1].match(/.{1,60}/g).join("\n")).join("\n")+"\n");
 fs.writeFileSync(`data/${name}_planted_truth.tsv`,"protein\tsignature\tstart\tend\n"+truth.join("\n")+"\n");
 let base=null;
 for(const [cn,ch,or] of CF)for(let rep=1;rep<=3;rep++){const r=run(prots,ch,or);const cov=new Set(r.all.map(l=>l.split("\t")[0])).size;
  bench.push([name,n,seed,cn,rep,r.ms.toFixed(1),r.rss.toFixed(1),r.all.length,(100*cov/n).toFixed(2),fp(r.all)].join(","));
  if(!base){base=r.all.slice().sort();fs.writeFileSync(`annotations/${name}_annotations.tsv`,base.join("\n")+"\n")}}
 const tset=new Set(truth.map(t=>{const f=t.split("\t");return f[0]+"\t"+f[1]}));const got=new Set(base.map(l=>{const f=l.split("\t");return f[0]+"\t"+f[4]}));
 let rec=0;tset.forEach(k=>{if(got.has(k))rec++});const cov=new Set(base.map(l=>l.split("\t")[0])).size;
 truthcmp.push([name,n,tset.size,rec,(100*rec/tset.size).toFixed(2),base.length,base.length-rec,(100*cov/n).toFixed(2)].join(","));
 for(const p of PAT){const h=base.filter(l=>l.split("\t")[4]===p[0]);bysig.push([name,p[0],'"'+p[1]+'"',h.length,new Set(h.map(l=>l.split("\t")[0])).size].join(","))}
 console.log("done",name)}
fs.writeFileSync("results/benchmark_runs.csv",bench.join("\n")+"\n");
fs.writeFileSync("results/planted_truth_comparison.csv",truthcmp.join("\n")+"\n");
fs.writeFileSync("results/matches_by_signature.csv",bysig.join("\n")+"\n");
const rows=bench.slice(1).map(l=>l.split(","));const sum=["dataset,runs,distinct_fingerprints,reproducible"];
for(const [n,seed] of SETS){const name=`synthetic_${n}_seed${seed}`,f=rows.filter(r=>r[0]===name).map(r=>r[9]);sum.push([name,f.length,new Set(f).size,new Set(f).size===1].join(","))}
fs.writeFileSync("results/reproducibility_summary.csv",sum.join("\n")+"\n");
