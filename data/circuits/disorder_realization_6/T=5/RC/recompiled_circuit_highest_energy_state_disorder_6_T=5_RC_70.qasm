OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1036086) q[0];
sx q[0];
rz(-1.684364) q[0];
sx q[0];
rz(1.3527704) q[0];
rz(-1.9569205) q[1];
sx q[1];
rz(-2.4689622) q[1];
sx q[1];
rz(-1.5157359) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11213451) q[0];
sx q[0];
rz(-1.3728276) q[0];
sx q[0];
rz(-2.9517118) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3358222) q[2];
sx q[2];
rz(-1.2122853) q[2];
sx q[2];
rz(1.1453978) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7616854) q[1];
sx q[1];
rz(-0.60098472) q[1];
sx q[1];
rz(-0.65324776) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8579086) q[3];
sx q[3];
rz(-1.9762044) q[3];
sx q[3];
rz(1.6943581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0036014) q[2];
sx q[2];
rz(-0.10819745) q[2];
sx q[2];
rz(-2.9555964) q[2];
rz(0.018639175) q[3];
sx q[3];
rz(-1.847307) q[3];
sx q[3];
rz(2.0712461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.5587392) q[0];
sx q[0];
rz(-2.5753729) q[0];
sx q[0];
rz(2.8232316) q[0];
rz(0.63703498) q[1];
sx q[1];
rz(-2.36167) q[1];
sx q[1];
rz(-0.46924082) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97142727) q[0];
sx q[0];
rz(-1.3999456) q[0];
sx q[0];
rz(2.9396482) q[0];
rz(-pi) q[1];
x q[1];
rz(0.083949371) q[2];
sx q[2];
rz(-1.9453887) q[2];
sx q[2];
rz(0.89886802) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7982895) q[1];
sx q[1];
rz(-1.8899916) q[1];
sx q[1];
rz(-1.809824) q[1];
rz(-pi) q[2];
rz(2.2607311) q[3];
sx q[3];
rz(-1.575205) q[3];
sx q[3];
rz(3.0772095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3176754) q[2];
sx q[2];
rz(-2.9313512) q[2];
sx q[2];
rz(2.4483185) q[2];
rz(-1.8215826) q[3];
sx q[3];
rz(-1.8157248) q[3];
sx q[3];
rz(2.0461953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0005118) q[0];
sx q[0];
rz(-0.63758272) q[0];
sx q[0];
rz(-0.47955036) q[0];
rz(0.50411049) q[1];
sx q[1];
rz(-1.9934374) q[1];
sx q[1];
rz(0.058301059) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8970015) q[0];
sx q[0];
rz(-1.0718597) q[0];
sx q[0];
rz(-1.8656436) q[0];
rz(-pi) q[1];
rz(-0.68966663) q[2];
sx q[2];
rz(-1.9885157) q[2];
sx q[2];
rz(0.4690241) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7229798) q[1];
sx q[1];
rz(-1.449391) q[1];
sx q[1];
rz(-1.0176119) q[1];
rz(-pi) q[2];
rz(0.4543484) q[3];
sx q[3];
rz(-0.55255167) q[3];
sx q[3];
rz(1.0187899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38640675) q[2];
sx q[2];
rz(-1.1815973) q[2];
sx q[2];
rz(0.03726658) q[2];
rz(-1.5197598) q[3];
sx q[3];
rz(-0.45791364) q[3];
sx q[3];
rz(0.38145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.195049) q[0];
sx q[0];
rz(-0.46849546) q[0];
sx q[0];
rz(-3.0750437) q[0];
rz(1.5244124) q[1];
sx q[1];
rz(-1.8981372) q[1];
sx q[1];
rz(0.7349416) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0733058) q[0];
sx q[0];
rz(-0.6319918) q[0];
sx q[0];
rz(0.39936622) q[0];
rz(-pi) q[1];
rz(-2.4114716) q[2];
sx q[2];
rz(-1.7459622) q[2];
sx q[2];
rz(-0.96792449) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0385602) q[1];
sx q[1];
rz(-1.8502697) q[1];
sx q[1];
rz(0.36906645) q[1];
rz(0.85020406) q[3];
sx q[3];
rz(-1.6805887) q[3];
sx q[3];
rz(-1.3444855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.28896004) q[2];
sx q[2];
rz(-1.5776878) q[2];
sx q[2];
rz(3.0328879) q[2];
rz(-0.92681256) q[3];
sx q[3];
rz(-0.97065297) q[3];
sx q[3];
rz(-1.577781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0573334) q[0];
sx q[0];
rz(-2.5045392) q[0];
sx q[0];
rz(-1.0908701) q[0];
rz(-0.57251656) q[1];
sx q[1];
rz(-1.9520091) q[1];
sx q[1];
rz(-2.3172839) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5474654) q[0];
sx q[0];
rz(-0.4484494) q[0];
sx q[0];
rz(-2.1350103) q[0];
rz(2.213843) q[2];
sx q[2];
rz(-1.130419) q[2];
sx q[2];
rz(-1.7945031) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5368176) q[1];
sx q[1];
rz(-1.0295492) q[1];
sx q[1];
rz(2.6763335) q[1];
rz(-pi) q[2];
rz(0.82309725) q[3];
sx q[3];
rz(-2.5783126) q[3];
sx q[3];
rz(-0.34264229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6606286) q[2];
sx q[2];
rz(-2.169256) q[2];
sx q[2];
rz(-0.29120293) q[2];
rz(-0.63699841) q[3];
sx q[3];
rz(-1.4642508) q[3];
sx q[3];
rz(1.8891107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3785192) q[0];
sx q[0];
rz(-2.851649) q[0];
sx q[0];
rz(-0.13105233) q[0];
rz(-1.1669) q[1];
sx q[1];
rz(-0.98375541) q[1];
sx q[1];
rz(-0.022620591) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2065434) q[0];
sx q[0];
rz(-0.66842043) q[0];
sx q[0];
rz(-2.5517625) q[0];
rz(-pi) q[1];
rz(-0.51949595) q[2];
sx q[2];
rz(-2.8970089) q[2];
sx q[2];
rz(-2.1313388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9359259) q[1];
sx q[1];
rz(-2.4956411) q[1];
sx q[1];
rz(-1.0420858) q[1];
rz(0.87225391) q[3];
sx q[3];
rz(-2.8917851) q[3];
sx q[3];
rz(1.2605309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.268198) q[2];
sx q[2];
rz(-2.4384273) q[2];
sx q[2];
rz(1.084729) q[2];
rz(-1.1449413) q[3];
sx q[3];
rz(-1.2866674) q[3];
sx q[3];
rz(1.0219319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5522083) q[0];
sx q[0];
rz(-1.5188058) q[0];
sx q[0];
rz(0.4735221) q[0];
rz(-1.9206958) q[1];
sx q[1];
rz(-2.7291606) q[1];
sx q[1];
rz(-1.268505) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6994239) q[0];
sx q[0];
rz(-2.0887362) q[0];
sx q[0];
rz(2.5479526) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67737214) q[2];
sx q[2];
rz(-2.2808135) q[2];
sx q[2];
rz(-1.7456919) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.699325) q[1];
sx q[1];
rz(-1.066477) q[1];
sx q[1];
rz(1.0902576) q[1];
rz(-pi) q[2];
rz(-0.082793391) q[3];
sx q[3];
rz(-2.0054545) q[3];
sx q[3];
rz(-2.5744048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1641757) q[2];
sx q[2];
rz(-1.5722534) q[2];
sx q[2];
rz(0.24641985) q[2];
rz(-0.94270802) q[3];
sx q[3];
rz(-1.0859414) q[3];
sx q[3];
rz(2.3781618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8134269) q[0];
sx q[0];
rz(-1.270371) q[0];
sx q[0];
rz(0.060591977) q[0];
rz(-1.6391485) q[1];
sx q[1];
rz(-1.7785347) q[1];
sx q[1];
rz(-1.5938119) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7029744) q[0];
sx q[0];
rz(-1.5697704) q[0];
sx q[0];
rz(0.0038504168) q[0];
rz(-pi) q[1];
rz(0.28569371) q[2];
sx q[2];
rz(-0.59425747) q[2];
sx q[2];
rz(0.93462925) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9309153) q[1];
sx q[1];
rz(-2.0077171) q[1];
sx q[1];
rz(-2.2411437) q[1];
rz(-pi) q[2];
rz(0.5129359) q[3];
sx q[3];
rz(-0.78425927) q[3];
sx q[3];
rz(-1.605214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.26555201) q[2];
sx q[2];
rz(-0.92090845) q[2];
sx q[2];
rz(1.457224) q[2];
rz(-0.11988104) q[3];
sx q[3];
rz(-0.20457743) q[3];
sx q[3];
rz(1.0286819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84751713) q[0];
sx q[0];
rz(-2.1471922) q[0];
sx q[0];
rz(1.1353528) q[0];
rz(-0.10336939) q[1];
sx q[1];
rz(-1.4048978) q[1];
sx q[1];
rz(0.9476544) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5948828) q[0];
sx q[0];
rz(-2.2077387) q[0];
sx q[0];
rz(1.4875436) q[0];
rz(2.060076) q[2];
sx q[2];
rz(-0.7502509) q[2];
sx q[2];
rz(2.3231151) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7458742) q[1];
sx q[1];
rz(-1.304879) q[1];
sx q[1];
rz(2.2661792) q[1];
rz(-pi) q[2];
rz(2.2470948) q[3];
sx q[3];
rz(-2.1966561) q[3];
sx q[3];
rz(1.8602399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.4244298) q[2];
sx q[2];
rz(-0.62817502) q[2];
sx q[2];
rz(-2.2372645) q[2];
rz(-1.2633911) q[3];
sx q[3];
rz(-1.2369913) q[3];
sx q[3];
rz(2.625107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.271027) q[0];
sx q[0];
rz(-2.3834383) q[0];
sx q[0];
rz(-0.98536056) q[0];
rz(2.789978) q[1];
sx q[1];
rz(-0.96013394) q[1];
sx q[1];
rz(-2.7947289) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064929334) q[0];
sx q[0];
rz(-2.4682625) q[0];
sx q[0];
rz(-1.2162186) q[0];
x q[1];
rz(3.0290481) q[2];
sx q[2];
rz(-2.0142609) q[2];
sx q[2];
rz(2.7751768) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.092791768) q[1];
sx q[1];
rz(-1.6597224) q[1];
sx q[1];
rz(2.4503178) q[1];
rz(-pi) q[2];
rz(1.204416) q[3];
sx q[3];
rz(-1.1561596) q[3];
sx q[3];
rz(2.7759027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.88006192) q[2];
sx q[2];
rz(-2.4226649) q[2];
sx q[2];
rz(3.0066709) q[2];
rz(-0.12923446) q[3];
sx q[3];
rz(-2.7415469) q[3];
sx q[3];
rz(-1.038704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1226817) q[0];
sx q[0];
rz(-1.2553348) q[0];
sx q[0];
rz(-3.0154764) q[0];
rz(0.46161721) q[1];
sx q[1];
rz(-1.4001662) q[1];
sx q[1];
rz(-2.9495159) q[1];
rz(0.89712894) q[2];
sx q[2];
rz(-1.9403602) q[2];
sx q[2];
rz(-0.42862949) q[2];
rz(0.57030525) q[3];
sx q[3];
rz(-2.8371596) q[3];
sx q[3];
rz(0.53573487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
