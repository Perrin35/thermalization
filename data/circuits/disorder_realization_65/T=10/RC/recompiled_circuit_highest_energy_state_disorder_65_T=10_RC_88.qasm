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
rz(0.49864545) q[0];
sx q[0];
rz(-2.5957624) q[0];
sx q[0];
rz(-0.11830615) q[0];
rz(-2.2220597) q[1];
sx q[1];
rz(-1.3032721) q[1];
sx q[1];
rz(-1.9272756) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34284232) q[0];
sx q[0];
rz(-1.3590706) q[0];
sx q[0];
rz(-1.1284853) q[0];
x q[1];
rz(2.5164175) q[2];
sx q[2];
rz(-2.0707651) q[2];
sx q[2];
rz(-1.2631877) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2976297) q[1];
sx q[1];
rz(-1.4342035) q[1];
sx q[1];
rz(0.94874391) q[1];
rz(-0.021539979) q[3];
sx q[3];
rz(-1.8439172) q[3];
sx q[3];
rz(2.5813951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.168657) q[2];
sx q[2];
rz(-2.0106222) q[2];
sx q[2];
rz(2.0449779) q[2];
rz(0.24122572) q[3];
sx q[3];
rz(-1.3696407) q[3];
sx q[3];
rz(-0.84478861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089652561) q[0];
sx q[0];
rz(-2.0515428) q[0];
sx q[0];
rz(2.7570778) q[0];
rz(-2.8334726) q[1];
sx q[1];
rz(-0.48079753) q[1];
sx q[1];
rz(-0.62517977) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2848672) q[0];
sx q[0];
rz(-1.856904) q[0];
sx q[0];
rz(1.5377783) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7169508) q[2];
sx q[2];
rz(-1.4745108) q[2];
sx q[2];
rz(1.7801628) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5601666) q[1];
sx q[1];
rz(-2.4871965) q[1];
sx q[1];
rz(-0.1398211) q[1];
rz(-pi) q[2];
rz(-2.1891648) q[3];
sx q[3];
rz(-1.980972) q[3];
sx q[3];
rz(0.71632121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1326617) q[2];
sx q[2];
rz(-0.93175685) q[2];
sx q[2];
rz(2.8399732) q[2];
rz(-0.53145069) q[3];
sx q[3];
rz(-0.88574946) q[3];
sx q[3];
rz(-1.1227932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6279491) q[0];
sx q[0];
rz(-1.0114089) q[0];
sx q[0];
rz(-1.1847786) q[0];
rz(0.14492497) q[1];
sx q[1];
rz(-1.8791608) q[1];
sx q[1];
rz(1.5941934) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7892706) q[0];
sx q[0];
rz(-1.0836067) q[0];
sx q[0];
rz(-2.4643023) q[0];
x q[1];
rz(-1.9504531) q[2];
sx q[2];
rz(-1.7413503) q[2];
sx q[2];
rz(1.7315239) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.84128296) q[1];
sx q[1];
rz(-2.1108859) q[1];
sx q[1];
rz(0.45068355) q[1];
rz(-1.2482485) q[3];
sx q[3];
rz(-1.1841906) q[3];
sx q[3];
rz(-2.4462388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3323815) q[2];
sx q[2];
rz(-2.1766162) q[2];
sx q[2];
rz(-1.3261718) q[2];
rz(2.2281846) q[3];
sx q[3];
rz(-1.923442) q[3];
sx q[3];
rz(-2.6083045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42295414) q[0];
sx q[0];
rz(-2.664743) q[0];
sx q[0];
rz(1.780321) q[0];
rz(-2.4286229) q[1];
sx q[1];
rz(-1.1402592) q[1];
sx q[1];
rz(0.68351173) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28302971) q[0];
sx q[0];
rz(-2.5682395) q[0];
sx q[0];
rz(1.5893776) q[0];
rz(1.7943125) q[2];
sx q[2];
rz(-1.7752675) q[2];
sx q[2];
rz(0.42753893) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5826368) q[1];
sx q[1];
rz(-0.92933944) q[1];
sx q[1];
rz(0.90822405) q[1];
x q[2];
rz(-2.8972647) q[3];
sx q[3];
rz(-1.4736946) q[3];
sx q[3];
rz(-1.3299143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5061364) q[2];
sx q[2];
rz(-0.41716245) q[2];
sx q[2];
rz(0.48640856) q[2];
rz(-2.6247315) q[3];
sx q[3];
rz(-1.9603399) q[3];
sx q[3];
rz(-2.1080871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085623398) q[0];
sx q[0];
rz(-1.6317246) q[0];
sx q[0];
rz(1.7686718) q[0];
rz(1.4068475) q[1];
sx q[1];
rz(-2.5739248) q[1];
sx q[1];
rz(0.47703201) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5236899) q[0];
sx q[0];
rz(-1.7672024) q[0];
sx q[0];
rz(1.1454885) q[0];
rz(-pi) q[1];
rz(-1.050694) q[2];
sx q[2];
rz(-1.2251266) q[2];
sx q[2];
rz(-3.1351948) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2998987) q[1];
sx q[1];
rz(-1.8704789) q[1];
sx q[1];
rz(1.0503286) q[1];
rz(-pi) q[2];
x q[2];
rz(2.414235) q[3];
sx q[3];
rz(-2.2148956) q[3];
sx q[3];
rz(1.3535045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.10739747) q[2];
sx q[2];
rz(-0.98661462) q[2];
sx q[2];
rz(2.6743215) q[2];
rz(-0.16921903) q[3];
sx q[3];
rz(-1.7305814) q[3];
sx q[3];
rz(-0.30538487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53448236) q[0];
sx q[0];
rz(-0.87962532) q[0];
sx q[0];
rz(2.95209) q[0];
rz(-2.8684008) q[1];
sx q[1];
rz(-2.0676282) q[1];
sx q[1];
rz(0.80064076) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9439745) q[0];
sx q[0];
rz(-1.6139784) q[0];
sx q[0];
rz(-1.3373242) q[0];
rz(-0.24836274) q[2];
sx q[2];
rz(-1.6585322) q[2];
sx q[2];
rz(1.1540292) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2792017) q[1];
sx q[1];
rz(-1.4420048) q[1];
sx q[1];
rz(-1.7679845) q[1];
rz(2.6012035) q[3];
sx q[3];
rz(-2.0107993) q[3];
sx q[3];
rz(-2.7203181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2022986) q[2];
sx q[2];
rz(-0.17387986) q[2];
sx q[2];
rz(-0.536971) q[2];
rz(-1.2823074) q[3];
sx q[3];
rz(-1.3064281) q[3];
sx q[3];
rz(-1.4793226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.110431) q[0];
sx q[0];
rz(-2.3261676) q[0];
sx q[0];
rz(-1.8102616) q[0];
rz(1.4415461) q[1];
sx q[1];
rz(-1.4794289) q[1];
sx q[1];
rz(1.2225245) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0095897) q[0];
sx q[0];
rz(-1.9870269) q[0];
sx q[0];
rz(0.25968174) q[0];
rz(-2.5897854) q[2];
sx q[2];
rz(-1.6344317) q[2];
sx q[2];
rz(-0.76630521) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4536087) q[1];
sx q[1];
rz(-2.2609432) q[1];
sx q[1];
rz(0.28120561) q[1];
rz(1.8172713) q[3];
sx q[3];
rz(-2.685084) q[3];
sx q[3];
rz(0.64303165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9407201) q[2];
sx q[2];
rz(-0.49786374) q[2];
sx q[2];
rz(-2.8087924) q[2];
rz(0.29916549) q[3];
sx q[3];
rz(-1.9009512) q[3];
sx q[3];
rz(-0.021473264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0126295) q[0];
sx q[0];
rz(-0.7681995) q[0];
sx q[0];
rz(0.27498883) q[0];
rz(3.0601652) q[1];
sx q[1];
rz(-0.80137253) q[1];
sx q[1];
rz(-1.3224695) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24182651) q[0];
sx q[0];
rz(-1.5112919) q[0];
sx q[0];
rz(-1.3391979) q[0];
rz(-pi) q[1];
rz(0.36825387) q[2];
sx q[2];
rz(-2.2241631) q[2];
sx q[2];
rz(-0.39619941) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6024633) q[1];
sx q[1];
rz(-2.2577598) q[1];
sx q[1];
rz(2.0571571) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4422936) q[3];
sx q[3];
rz(-0.85996503) q[3];
sx q[3];
rz(3.1131203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4203809) q[2];
sx q[2];
rz(-0.12004852) q[2];
sx q[2];
rz(-1.7397286) q[2];
rz(-1.8574572) q[3];
sx q[3];
rz(-1.9522342) q[3];
sx q[3];
rz(0.0349667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5597124) q[0];
sx q[0];
rz(-1.548883) q[0];
sx q[0];
rz(2.6522719) q[0];
rz(3.1210461) q[1];
sx q[1];
rz(-1.2362044) q[1];
sx q[1];
rz(0.70125854) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1398292) q[0];
sx q[0];
rz(-2.1368623) q[0];
sx q[0];
rz(2.8553814) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0938265) q[2];
sx q[2];
rz(-0.86776185) q[2];
sx q[2];
rz(0.060279708) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32966954) q[1];
sx q[1];
rz(-0.81336248) q[1];
sx q[1];
rz(1.2438891) q[1];
rz(-pi) q[2];
rz(-0.55072983) q[3];
sx q[3];
rz(-0.76590532) q[3];
sx q[3];
rz(0.85635161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7265085) q[2];
sx q[2];
rz(-0.61843151) q[2];
sx q[2];
rz(-0.48039594) q[2];
rz(-1.4091617) q[3];
sx q[3];
rz(-1.3536072) q[3];
sx q[3];
rz(-2.8044146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3163863) q[0];
sx q[0];
rz(-1.9012863) q[0];
sx q[0];
rz(-0.41148841) q[0];
rz(-1.8972137) q[1];
sx q[1];
rz(-2.0020516) q[1];
sx q[1];
rz(0.32434514) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2712944) q[0];
sx q[0];
rz(-0.37640262) q[0];
sx q[0];
rz(-1.7145304) q[0];
x q[1];
rz(0.34180832) q[2];
sx q[2];
rz(-1.2380307) q[2];
sx q[2];
rz(1.8012799) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.577475) q[1];
sx q[1];
rz(-1.0077969) q[1];
sx q[1];
rz(-0.098280829) q[1];
rz(-1.5839229) q[3];
sx q[3];
rz(-1.6919976) q[3];
sx q[3];
rz(1.2016202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1856498) q[2];
sx q[2];
rz(-2.0273384) q[2];
sx q[2];
rz(-0.11697098) q[2];
rz(-2.1864435) q[3];
sx q[3];
rz(-2.9626466) q[3];
sx q[3];
rz(-2.5999033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2627926) q[0];
sx q[0];
rz(-2.265082) q[0];
sx q[0];
rz(-2.2646917) q[0];
rz(1.1204489) q[1];
sx q[1];
rz(-0.81605492) q[1];
sx q[1];
rz(1.1647404) q[1];
rz(1.9399613) q[2];
sx q[2];
rz(-2.6045447) q[2];
sx q[2];
rz(2.4640026) q[2];
rz(-1.8414303) q[3];
sx q[3];
rz(-1.5951372) q[3];
sx q[3];
rz(-1.5230877) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
