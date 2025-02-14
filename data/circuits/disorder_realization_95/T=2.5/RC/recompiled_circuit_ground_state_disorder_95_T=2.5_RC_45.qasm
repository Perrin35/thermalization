OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9606544) q[0];
sx q[0];
rz(-0.063194312) q[0];
sx q[0];
rz(-0.59046459) q[0];
rz(-0.2904627) q[1];
sx q[1];
rz(1.7658748) q[1];
sx q[1];
rz(9.4212846) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1336023) q[0];
sx q[0];
rz(-0.65235814) q[0];
sx q[0];
rz(-2.2732387) q[0];
rz(2.2512746) q[2];
sx q[2];
rz(-1.3290452) q[2];
sx q[2];
rz(1.7861451) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0776135) q[1];
sx q[1];
rz(-3.1099456) q[1];
sx q[1];
rz(2.3961303) q[1];
rz(-2.6925647) q[3];
sx q[3];
rz(-1.0211583) q[3];
sx q[3];
rz(0.47269905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8840088) q[2];
sx q[2];
rz(-0.85870063) q[2];
sx q[2];
rz(-1.0165455) q[2];
rz(-2.2782585) q[3];
sx q[3];
rz(-0.038067929) q[3];
sx q[3];
rz(-1.6421002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4676374) q[0];
sx q[0];
rz(-1.4731982) q[0];
sx q[0];
rz(3.0067645) q[0];
rz(-0.22678953) q[1];
sx q[1];
rz(-0.062954523) q[1];
sx q[1];
rz(-0.24980587) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733179) q[0];
sx q[0];
rz(-2.3767243) q[0];
sx q[0];
rz(-2.5165783) q[0];
rz(-pi) q[1];
x q[1];
rz(0.05257865) q[2];
sx q[2];
rz(-2.0705531) q[2];
sx q[2];
rz(1.6249958) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0909101) q[1];
sx q[1];
rz(-1.5600648) q[1];
sx q[1];
rz(-0.077535943) q[1];
x q[2];
rz(-3.0071179) q[3];
sx q[3];
rz(-1.2081097) q[3];
sx q[3];
rz(0.69512007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2251542) q[2];
sx q[2];
rz(-0.068000451) q[2];
sx q[2];
rz(0.98006836) q[2];
rz(0.89503908) q[3];
sx q[3];
rz(-0.74256247) q[3];
sx q[3];
rz(1.1802973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5880244) q[0];
sx q[0];
rz(-0.50694412) q[0];
sx q[0];
rz(-1.5356327) q[0];
rz(-1.5060679) q[1];
sx q[1];
rz(-0.82553828) q[1];
sx q[1];
rz(0.95605409) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0900537) q[0];
sx q[0];
rz(-2.9057626) q[0];
sx q[0];
rz(2.7084742) q[0];
rz(-pi) q[1];
rz(1.6371085) q[2];
sx q[2];
rz(-2.300718) q[2];
sx q[2];
rz(-0.17375565) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.07900958) q[1];
sx q[1];
rz(-2.2979567) q[1];
sx q[1];
rz(0.14240245) q[1];
rz(-pi) q[2];
rz(-2.6272229) q[3];
sx q[3];
rz(-1.2651431) q[3];
sx q[3];
rz(2.5955868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8354336) q[2];
sx q[2];
rz(-0.71194887) q[2];
sx q[2];
rz(-0.63208675) q[2];
rz(0.88405526) q[3];
sx q[3];
rz(-0.017280936) q[3];
sx q[3];
rz(-2.250905) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66568351) q[0];
sx q[0];
rz(-1.8432239) q[0];
sx q[0];
rz(-0.16715288) q[0];
rz(-1.2627603) q[1];
sx q[1];
rz(-2.1985168) q[1];
sx q[1];
rz(-1.7335588) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4830925) q[0];
sx q[0];
rz(-1.8611055) q[0];
sx q[0];
rz(2.6997552) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7825259) q[2];
sx q[2];
rz(-1.5671726) q[2];
sx q[2];
rz(-0.45932367) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7552897) q[1];
sx q[1];
rz(-1.7610671) q[1];
sx q[1];
rz(-2.931837) q[1];
x q[2];
rz(-1.2795598) q[3];
sx q[3];
rz(-2.3197745) q[3];
sx q[3];
rz(-2.237221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.35272804) q[2];
sx q[2];
rz(-0.015268607) q[2];
sx q[2];
rz(-0.35066476) q[2];
rz(0.26385012) q[3];
sx q[3];
rz(-0.00084547384) q[3];
sx q[3];
rz(-1.2034169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26461399) q[0];
sx q[0];
rz(-2.3961841) q[0];
sx q[0];
rz(1.4234446) q[0];
rz(0.28165948) q[1];
sx q[1];
rz(-1.6335082) q[1];
sx q[1];
rz(-1.8219832) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9284009) q[0];
sx q[0];
rz(-1.9921148) q[0];
sx q[0];
rz(-2.1405959) q[0];
x q[1];
rz(-1.5983888) q[2];
sx q[2];
rz(-1.5910307) q[2];
sx q[2];
rz(-1.2638365) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4180243) q[1];
sx q[1];
rz(-0.71360525) q[1];
sx q[1];
rz(0.51939555) q[1];
rz(1.237147) q[3];
sx q[3];
rz(-1.1859484) q[3];
sx q[3];
rz(2.9001482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.264297) q[2];
sx q[2];
rz(-3.0946315) q[2];
sx q[2];
rz(-1.7128672) q[2];
rz(1.9127539) q[3];
sx q[3];
rz(-0.052534025) q[3];
sx q[3];
rz(0.088168941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0577211) q[0];
sx q[0];
rz(-0.11744048) q[0];
sx q[0];
rz(-0.55026662) q[0];
rz(-0.30217198) q[1];
sx q[1];
rz(-1.7487532) q[1];
sx q[1];
rz(-2.7064533) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36212039) q[0];
sx q[0];
rz(-1.605662) q[0];
sx q[0];
rz(-2.930234) q[0];
x q[1];
rz(-3.1341564) q[2];
sx q[2];
rz(-1.5698293) q[2];
sx q[2];
rz(-1.8966248) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8813215) q[1];
sx q[1];
rz(-1.146293) q[1];
sx q[1];
rz(-3.0851302) q[1];
x q[2];
rz(-2.4184259) q[3];
sx q[3];
rz(-2.1373014) q[3];
sx q[3];
rz(0.97149937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90007323) q[2];
sx q[2];
rz(-3.1406431) q[2];
sx q[2];
rz(2.1692236) q[2];
rz(2.1420245) q[3];
sx q[3];
rz(-3.1344423) q[3];
sx q[3];
rz(2.0991367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5480176) q[0];
sx q[0];
rz(-1.2418208) q[0];
sx q[0];
rz(-1.0663363) q[0];
rz(-1.713133) q[1];
sx q[1];
rz(-2.8722873) q[1];
sx q[1];
rz(1.8535293) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7748874) q[0];
sx q[0];
rz(-0.88434356) q[0];
sx q[0];
rz(1.9839601) q[0];
rz(-pi) q[1];
rz(0.84128235) q[2];
sx q[2];
rz(-1.3026574) q[2];
sx q[2];
rz(-0.27716217) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8180698) q[1];
sx q[1];
rz(-1.5574291) q[1];
sx q[1];
rz(-1.6105777) q[1];
rz(-pi) q[2];
rz(-2.101048) q[3];
sx q[3];
rz(-1.7826579) q[3];
sx q[3];
rz(-2.7303641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.11918934) q[2];
sx q[2];
rz(-0.069312118) q[2];
sx q[2];
rz(2.8663087) q[2];
rz(-2.9149808) q[3];
sx q[3];
rz(-0.35609069) q[3];
sx q[3];
rz(0.4376469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3771123) q[0];
sx q[0];
rz(-2.8545916) q[0];
sx q[0];
rz(2.5544033) q[0];
rz(1.6181234) q[1];
sx q[1];
rz(-1.1174997) q[1];
sx q[1];
rz(-1.5709741) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.068741) q[0];
sx q[0];
rz(-0.91813164) q[0];
sx q[0];
rz(1.0314343) q[0];
rz(1.1067939) q[2];
sx q[2];
rz(-1.6308068) q[2];
sx q[2];
rz(-1.6911708) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5779004) q[1];
sx q[1];
rz(-1.6017388) q[1];
sx q[1];
rz(-1.5269431) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0944902) q[3];
sx q[3];
rz(-2.5137551) q[3];
sx q[3];
rz(-1.353457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9470584) q[2];
sx q[2];
rz(-2.9583866) q[2];
sx q[2];
rz(1.3018695) q[2];
rz(0.434508) q[3];
sx q[3];
rz(-0.040381581) q[3];
sx q[3];
rz(0.44809189) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3473564) q[0];
sx q[0];
rz(-0.11581049) q[0];
sx q[0];
rz(-1.7606803) q[0];
rz(1.5877089) q[1];
sx q[1];
rz(-0.96276182) q[1];
sx q[1];
rz(0.10996058) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9048683) q[0];
sx q[0];
rz(-2.3873608) q[0];
sx q[0];
rz(-1.8484902) q[0];
rz(-pi) q[1];
rz(-0.25935632) q[2];
sx q[2];
rz(-2.4217941) q[2];
sx q[2];
rz(0.46611818) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.0011855652) q[1];
sx q[1];
rz(-2.3275536) q[1];
sx q[1];
rz(-3.1073529) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4110097) q[3];
sx q[3];
rz(-1.8786123) q[3];
sx q[3];
rz(-2.6871443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.79003698) q[2];
sx q[2];
rz(-2.07708) q[2];
sx q[2];
rz(2.7891187) q[2];
rz(0.6414837) q[3];
sx q[3];
rz(-3.0951169) q[3];
sx q[3];
rz(-0.36267734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27956692) q[0];
sx q[0];
rz(-0.27722219) q[0];
sx q[0];
rz(-0.56420457) q[0];
rz(1.5203681) q[1];
sx q[1];
rz(-1.0313326) q[1];
sx q[1];
rz(0.079781562) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4046116) q[0];
sx q[0];
rz(-0.68697646) q[0];
sx q[0];
rz(-1.8115329) q[0];
rz(-1.6185845) q[2];
sx q[2];
rz(-1.6182247) q[2];
sx q[2];
rz(-0.090605926) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33377117) q[1];
sx q[1];
rz(-0.83955315) q[1];
sx q[1];
rz(0.34656406) q[1];
rz(2.8989927) q[3];
sx q[3];
rz(-1.4897457) q[3];
sx q[3];
rz(-1.9417861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.12678777) q[2];
sx q[2];
rz(-0.0099651907) q[2];
sx q[2];
rz(1.0753746) q[2];
rz(2.3672095) q[3];
sx q[3];
rz(-0.024024809) q[3];
sx q[3];
rz(-0.27124852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8423691) q[0];
sx q[0];
rz(-1.5726226) q[0];
sx q[0];
rz(1.5725305) q[0];
rz(-2.6161999) q[1];
sx q[1];
rz(-3.0699984) q[1];
sx q[1];
rz(2.933554) q[1];
rz(-2.9276949) q[2];
sx q[2];
rz(-2.1281617) q[2];
sx q[2];
rz(-2.8585363) q[2];
rz(2.0490859) q[3];
sx q[3];
rz(-2.1989965) q[3];
sx q[3];
rz(0.60529706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
