// Real figures transcribed from Blum et al. 2026, Bioinformatics Advances 6(1):vbag141 (CC BY 4.0), Tables 1 and 2 and Section 3.2.
const fs=require('fs');
const T1=[["A. thaliana","UP000006548",27448,"9 h 29 min",569,"6 h 10 min",370,5],["C. elegans","UP000001940",19824,"5 h 30 min",330,"3 h 49 min",229,6],["D. melanogaster","UP000000803",13822,"5 h 59 min",359,"4 h 14 min",254,6],["E. coli","UP000000625",4402,"1 h 29 min",89,"0 h 39 min",39,2],["H. sapiens","UP000005640",20659,"12 h 13 min",733,"6 h 42 min",402,14],["M. musculus","UP000000589",21803,"11 h 33 min",693,"6 h 41 min",401,13],["M. tuberculosis","UP000001584",3996,"1 h 39 min",99,"0 h 44 min",44,2],["S. cerevisiae","UP000002311",6065,"2 h 54 min",174,"1 h 22 min",82,3],["Z. mays","UP000007305",39209,"10 h 51 min",651,"7 h 28 min",448,5]];
const T2=[["CDD",1,1,1],["CATH-Gene3D",.998,.999,.998],["HAMAP",1,1,1],["NCBIFAM",1,1,1],["PANTHER",.985,.999,.992],["Pfam",1,1,1],["PIRSF",.999,.999,.999],["PRINTS",1,1,1],["PROSITE Patterns",.999,.999,.999],["PROSITE Profiles",1,1,1],["SFLD",1,1,1],["SMART",.999,.993,.996],["SUPERFAMILY",1,1,1]];
const CL=[["AWS","single instance",10.14],["Google Cloud","single instance",8.13],["AWS","hybrid (Batch, spot)",12.43],["Google Cloud","hybrid (Batch, spot)",11.53]];
fs.writeFileSync("results/paper_table1_runtimes.csv","proteome,uniprot_id,sequences,ips5_reported,ips5_min,ips6_reported,ips6_min,ips6_plus_match_lookup_min,ips6_reduction_pct_computed,speedup_x_computed\n"+T1.map(r=>[r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],(100*(1-r[6]/r[4])).toFixed(1),(r[4]/r[6]).toFixed(2)].join(",")).join("\n")+"\n");
fs.writeFileSync("results/paper_table2_concordance.csv","member_database,precision,sensitivity,f1\n"+T2.map(r=>r.join(",")).join("\n")+"\n");
fs.writeFileSync("results/paper_cloud_costs_human_proteome.csv","provider,mode,cost_usd\n"+CL.map(r=>r.join(",")).join("\n")+"\n");
const max=733,bar=(v,c)=>`<span style="display:block;height:9px;border-radius:2px;background:var(${c});width:${(100*v/max).toFixed(1)}%;margin:2px 0"></span>`;
const t1=`<div class="wrap"><table><tr><th>Proteome</th><th>Sequences</th><th>IPS 5</th><th>IPS 6</th><th>IPS 6 + match lookup</th><th>IPS 6 time saved*</th><th style="min-width:160px">IPS 5 (orange) vs IPS 6 (teal)</th></tr>`+T1.map(r=>`<tr><td><i>${r[0]}</i><br><small>${r[1]}</small></td><td>${r[2].toLocaleString("en-GB")}</td><td>${r[3]}</td><td>${r[5]}</td><td>${r[7]} min</td><td>${(100*(1-r[6]/r[4])).toFixed(0)}% (${(r[4]/r[6]).toFixed(1)}x)</td><td>${bar(r[4],"--b")}${bar(r[6],"--a")}</td></tr>`).join("")+`</table></div><p class="sub" style="font-size:.9rem">*Computed here from the reported times, not stated in the paper in this form.</p>`;
const t2=`<div class="wrap"><table><tr><th>Member database</th><th>Precision</th><th>Sensitivity</th><th>F1</th></tr>`+T2.map(r=>`<tr><td>${r[0]}</td><td>${r[1]}</td><td>${r[2]}</td><td>${r[3]}</td></tr>`).join("")+`</table></div>`;
const t3=`<div class="wrap"><table><tr><th>Provider</th><th>Mode</th><th>Cost, human proteome</th></tr>`+CL.map(r=>`<tr><td>${r[0]}</td><td>${r[1]}</td><td>$${r[2].toFixed(2)}</td></tr>`).join("")+`</table></div>`;
const sec=`<section data-p="/paper" hidden>
<h1>Real data from the paper</h1>
<p>These are the published results of InterProScan 6, copied from Blum M, Hobbs E, Florentino L, Bateman A. <i>InterProScan 6: a modern large-scale protein function annotation pipeline.</i> Bioinformatics Advances 2026; 6(1): vbag141. <a href="https://doi.org/10.1093/bioadv/vbag141" style="color:var(--a)">doi:10.1093/bioadv/vbag141</a> (open access, CC BY 4.0). I did not run anything to get them.</p>
<h2>How the paper measured it</h2>
<ul><li>Nine UniProt reference proteomes (release 2025_04, one protein sequence per gene).</li>
<li>InterProScan 5.77-108.0 against InterProScan 6.0.0 run through Nextflow with Singularity, both using InterPro release 108.0.</li>
<li>EMBL-EBI Slurm cluster, 8 CPU cores per task (Intel Xeon Gold 6252), 90 GB memory requested per job.</li>
<li>Median of five runs. Times are full pipeline runs and exclude one-time data downloads. Run-to-run spread (MAD/median) was under 6%.</li></ul>
<h2>Runtime, Table 1</h2>
${t1}
<p>With the match lookup switched on, sequences already known to UniParc are looked up by checksum instead of being analysed again, which is why those times are minutes.</p>
<h2>Agreement between versions, Table 2</h2>
<p>All 573,661 Swiss-Prot sequences (release 2025_04). InterProScan 5 is the reference, and two matches count as the same when they overlap reciprocally by at least 80%.</p>
${t2}
<h2>Cloud cost, Section 3.2</h2>
<p>Full human proteome. The hybrid runs took under two hours and used spot instances.</p>
${t3}
<h2>Using this for your reproduction</h2>
<ul><li>The two smallest proteomes are the realistic ones to repeat: <i>E. coli</i> (UP000000625, 4,402 sequences, 39 min in the paper) and <i>M. tuberculosis</i> (UP000001584, 3,996 sequences, 44 min). Download them from UniProt using the exact IDs above.</li>
<li>Your hardware will differ, so compare ratios (IPS 6 with and without match lookup) rather than raw minutes, and state your CPU count.</li>
<li>Use the same InterPro release as the paper (108.0) when checking annotations, so database changes do not confuse the comparison.</li>
<li>The paper also lists every individual benchmark run in its Supplementary Table 1 (an xlsx file on the journal page), which I could not download.</li></ul>
<p>Machine-readable copies of these tables are in the zip under <code>results/paper_*.csv</code>.</p>
</section>

`;
fs.writeFileSync("/tmp/paper_section.html",sec);
