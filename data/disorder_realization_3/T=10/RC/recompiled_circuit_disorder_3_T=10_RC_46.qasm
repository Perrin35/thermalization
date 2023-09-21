OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(2.6565235) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(0.41419849) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6611377) q[0];
sx q[0];
rz(-1.9248795) q[0];
sx q[0];
rz(2.7829091) q[0];
rz(-pi) q[1];
rz(-3.0468416) q[2];
sx q[2];
rz(-1.0031327) q[2];
sx q[2];
rz(1.809158) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58798446) q[1];
sx q[1];
rz(-1.6349287) q[1];
sx q[1];
rz(1.1448121) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26478404) q[3];
sx q[3];
rz(-1.6695108) q[3];
sx q[3];
rz(-0.36987723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2177314) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(-0.031575354) q[2];
rz(-1.8850373) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(-0.1698499) q[0];
rz(2.4376712) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(2.6020715) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8115494) q[0];
sx q[0];
rz(-2.0199611) q[0];
sx q[0];
rz(-3.0898068) q[0];
rz(2.6846634) q[2];
sx q[2];
rz(-2.3655342) q[2];
sx q[2];
rz(-1.4571783) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62944618) q[1];
sx q[1];
rz(-2.0680032) q[1];
sx q[1];
rz(0.21558233) q[1];
x q[2];
rz(-0.94445618) q[3];
sx q[3];
rz(-0.8144905) q[3];
sx q[3];
rz(-2.8163547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4743621) q[2];
sx q[2];
rz(-1.903406) q[2];
sx q[2];
rz(-0.24307069) q[2];
rz(2.4754751) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(-1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77984017) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(-0.4483805) q[0];
rz(1.386863) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(-0.2562491) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1218131) q[0];
sx q[0];
rz(-2.461713) q[0];
sx q[0];
rz(-2.7133184) q[0];
x q[1];
rz(-2.3149812) q[2];
sx q[2];
rz(-0.6140784) q[2];
sx q[2];
rz(2.5342864) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5445404) q[1];
sx q[1];
rz(-1.8982732) q[1];
sx q[1];
rz(-2.271133) q[1];
rz(-pi) q[2];
rz(0.17685299) q[3];
sx q[3];
rz(-1.5336509) q[3];
sx q[3];
rz(-1.1249441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8024575) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(0.57470542) q[2];
rz(-1.7859219) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(-1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9280076) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(-2.8821049) q[0];
rz(-1.150594) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(-2.4096699) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61011945) q[0];
sx q[0];
rz(-2.0391132) q[0];
sx q[0];
rz(-2.5584695) q[0];
rz(-1.6875793) q[2];
sx q[2];
rz(-2.2846966) q[2];
sx q[2];
rz(1.2920979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.96958292) q[1];
sx q[1];
rz(-1.828555) q[1];
sx q[1];
rz(1.7715363) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6287574) q[3];
sx q[3];
rz(-2.8898015) q[3];
sx q[3];
rz(-0.73392111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6966454) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(0.15110061) q[2];
rz(-2.5949196) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(-2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5916409) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(1.8792101) q[0];
rz(1.6732015) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(-0.79777065) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3691468) q[0];
sx q[0];
rz(-1.4506842) q[0];
sx q[0];
rz(3.1390879) q[0];
rz(-pi) q[1];
rz(-1.351379) q[2];
sx q[2];
rz(-1.8674388) q[2];
sx q[2];
rz(-0.23362939) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5835727) q[1];
sx q[1];
rz(-2.1961374) q[1];
sx q[1];
rz(-3.0261092) q[1];
x q[2];
rz(0.62545125) q[3];
sx q[3];
rz(-2.6126385) q[3];
sx q[3];
rz(1.7367425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.795934) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(-2.8395555) q[2];
rz(-1.1473514) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(-2.5938477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40925947) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(1.1557895) q[0];
rz(-1.0844768) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(0.070080431) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4159211) q[0];
sx q[0];
rz(-1.3999181) q[0];
sx q[0];
rz(0.15539774) q[0];
x q[1];
rz(-0.55307936) q[2];
sx q[2];
rz(-2.2773909) q[2];
sx q[2];
rz(-1.7318219) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.97092123) q[1];
sx q[1];
rz(-1.0075924) q[1];
sx q[1];
rz(-1.3794273) q[1];
x q[2];
rz(-1.1280941) q[3];
sx q[3];
rz(-0.08249313) q[3];
sx q[3];
rz(-3.1174297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8391116) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(-2.5202259) q[2];
rz(1.4403884) q[3];
sx q[3];
rz(-0.50783235) q[3];
sx q[3];
rz(2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086534111) q[0];
sx q[0];
rz(-1.7250412) q[0];
sx q[0];
rz(-2.281718) q[0];
rz(1.2043918) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(-0.0079356114) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7746349) q[0];
sx q[0];
rz(-1.3152221) q[0];
sx q[0];
rz(0.5704244) q[0];
rz(-pi) q[1];
rz(-2.2631049) q[2];
sx q[2];
rz(-1.2307067) q[2];
sx q[2];
rz(-2.0932587) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4850033) q[1];
sx q[1];
rz(-2.2409229) q[1];
sx q[1];
rz(-0.90799241) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8273724) q[3];
sx q[3];
rz(-1.946297) q[3];
sx q[3];
rz(0.34208959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4456711) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(2.725214) q[2];
rz(1.773206) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96173441) q[0];
sx q[0];
rz(-0.27357736) q[0];
sx q[0];
rz(-0.36488786) q[0];
rz(0.94003135) q[1];
sx q[1];
rz(-2.6049728) q[1];
sx q[1];
rz(-1.5023124) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84425981) q[0];
sx q[0];
rz(-1.5607921) q[0];
sx q[0];
rz(-2.8914333) q[0];
rz(2.3055326) q[2];
sx q[2];
rz(-1.9327455) q[2];
sx q[2];
rz(-1.8723633) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.18174905) q[1];
sx q[1];
rz(-2.1310924) q[1];
sx q[1];
rz(-0.77397857) q[1];
rz(-pi) q[2];
rz(-1.3650465) q[3];
sx q[3];
rz(-1.660941) q[3];
sx q[3];
rz(2.3004325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49729785) q[2];
sx q[2];
rz(-0.50540322) q[2];
sx q[2];
rz(2.0765182) q[2];
rz(2.8403357) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(-1.5208972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168468) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(-2.705943) q[0];
rz(1.3849974) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(-0.41697821) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5676253) q[0];
sx q[0];
rz(-1.6084533) q[0];
sx q[0];
rz(0.91659878) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4776393) q[2];
sx q[2];
rz(-2.0400527) q[2];
sx q[2];
rz(2.9479153) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4638085) q[1];
sx q[1];
rz(-0.52075547) q[1];
sx q[1];
rz(-1.2893454) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0005089) q[3];
sx q[3];
rz(-1.699563) q[3];
sx q[3];
rz(1.5537378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0372662) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(-0.22582516) q[2];
rz(2.9337692) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(-2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.726783) q[0];
sx q[0];
rz(-0.2668969) q[0];
sx q[0];
rz(-1.5243994) q[0];
rz(-0.95364755) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(-1.3226002) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029862558) q[0];
sx q[0];
rz(-0.39806453) q[0];
sx q[0];
rz(-1.2533305) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0820504) q[2];
sx q[2];
rz(-2.6419123) q[2];
sx q[2];
rz(-0.79364712) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9433371) q[1];
sx q[1];
rz(-0.88611929) q[1];
sx q[1];
rz(-0.68590045) q[1];
rz(-0.91536509) q[3];
sx q[3];
rz(-1.4487106) q[3];
sx q[3];
rz(-1.9742427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7913197) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(-2.1255169) q[2];
rz(-1.919205) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(-0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14324698) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(1.3636419) q[1];
sx q[1];
rz(-1.9121871) q[1];
sx q[1];
rz(-1.8180064) q[1];
rz(-2.6323742) q[2];
sx q[2];
rz(-1.5794532) q[2];
sx q[2];
rz(3.0184359) q[2];
rz(-1.3397459) q[3];
sx q[3];
rz(-1.6587202) q[3];
sx q[3];
rz(-0.29826577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];