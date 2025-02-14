OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8483228) q[0];
sx q[0];
rz(-0.25265101) q[0];
sx q[0];
rz(1.9360315) q[0];
rz(-1.128101) q[1];
sx q[1];
rz(4.4681273) q[1];
sx q[1];
rz(11.820643) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099961258) q[0];
sx q[0];
rz(-1.7929165) q[0];
sx q[0];
rz(-3.0801386) q[0];
rz(-0.029736515) q[2];
sx q[2];
rz(-2.9171037) q[2];
sx q[2];
rz(-2.9475074) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.28994014) q[1];
sx q[1];
rz(-1.9059685) q[1];
sx q[1];
rz(0.22214684) q[1];
rz(-1.0966572) q[3];
sx q[3];
rz(-1.7563987) q[3];
sx q[3];
rz(1.4193648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2241609) q[2];
sx q[2];
rz(-1.5956343) q[2];
sx q[2];
rz(1.4707627) q[2];
rz(-1.6338232) q[3];
sx q[3];
rz(-0.016851146) q[3];
sx q[3];
rz(2.2182218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725752) q[0];
sx q[0];
rz(-1.9439531) q[0];
sx q[0];
rz(-1.5684599) q[0];
rz(-2.9720427) q[1];
sx q[1];
rz(-0.1145656) q[1];
sx q[1];
rz(0.13470185) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54808211) q[0];
sx q[0];
rz(-2.5786886) q[0];
sx q[0];
rz(-0.68771945) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1293389) q[2];
sx q[2];
rz(-1.6394221) q[2];
sx q[2];
rz(1.5315646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4410494) q[1];
sx q[1];
rz(-1.8722539) q[1];
sx q[1];
rz(-0.45977199) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.080218519) q[3];
sx q[3];
rz(-1.3933946) q[3];
sx q[3];
rz(-2.3631848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0160825) q[2];
sx q[2];
rz(-1.4849911) q[2];
sx q[2];
rz(2.9888195) q[2];
rz(1.3656535) q[3];
sx q[3];
rz(-3.1047265) q[3];
sx q[3];
rz(-0.16807817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50853866) q[0];
sx q[0];
rz(-2.3336053) q[0];
sx q[0];
rz(-0.48164865) q[0];
rz(2.956849) q[1];
sx q[1];
rz(-1.7748723) q[1];
sx q[1];
rz(0.99536037) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48993123) q[0];
sx q[0];
rz(-1.7999987) q[0];
sx q[0];
rz(1.9263173) q[0];
rz(0.048599343) q[2];
sx q[2];
rz(-1.6311797) q[2];
sx q[2];
rz(0.27602613) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9650594) q[1];
sx q[1];
rz(-0.73672026) q[1];
sx q[1];
rz(0.45335575) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97484346) q[3];
sx q[3];
rz(-1.5454834) q[3];
sx q[3];
rz(1.5060177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.92622009) q[2];
sx q[2];
rz(-0.054457713) q[2];
sx q[2];
rz(-0.064662956) q[2];
rz(1.0489382) q[3];
sx q[3];
rz(-3.1147396) q[3];
sx q[3];
rz(1.2803199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8365086) q[0];
sx q[0];
rz(-3.0236112) q[0];
sx q[0];
rz(-0.82103658) q[0];
rz(3.0631284) q[1];
sx q[1];
rz(-1.4469701) q[1];
sx q[1];
rz(0.99748126) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067046384) q[0];
sx q[0];
rz(-2.1830478) q[0];
sx q[0];
rz(1.1825829) q[0];
rz(-pi) q[1];
rz(3.1053379) q[2];
sx q[2];
rz(-1.5116351) q[2];
sx q[2];
rz(-1.349337) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9033244) q[1];
sx q[1];
rz(-2.0966171) q[1];
sx q[1];
rz(-2.8366798) q[1];
rz(-pi) q[2];
rz(2.663199) q[3];
sx q[3];
rz(-1.5198623) q[3];
sx q[3];
rz(-0.8138322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0492101) q[2];
sx q[2];
rz(-3.1123078) q[2];
sx q[2];
rz(-1.9878261) q[2];
rz(0.18618259) q[3];
sx q[3];
rz(-3.0598873) q[3];
sx q[3];
rz(0.33128273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1368197) q[0];
sx q[0];
rz(-0.81757075) q[0];
sx q[0];
rz(-1.9983043) q[0];
rz(1.1413057) q[1];
sx q[1];
rz(-0.80991304) q[1];
sx q[1];
rz(0.51796651) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2835613) q[0];
sx q[0];
rz(-1.0378583) q[0];
sx q[0];
rz(0.74646797) q[0];
x q[1];
rz(1.632156) q[2];
sx q[2];
rz(-0.01505919) q[2];
sx q[2];
rz(1.3850152) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.28455602) q[1];
sx q[1];
rz(-1.9003881) q[1];
sx q[1];
rz(0.28917851) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44995309) q[3];
sx q[3];
rz(-0.3287238) q[3];
sx q[3];
rz(2.4462819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3190069) q[2];
sx q[2];
rz(-2.1876882) q[2];
sx q[2];
rz(0.32773584) q[2];
rz(1.94708) q[3];
sx q[3];
rz(-0.12735282) q[3];
sx q[3];
rz(-2.2849042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0026534) q[0];
sx q[0];
rz(-2.8728573) q[0];
sx q[0];
rz(-0.56185454) q[0];
rz(1.4909164) q[1];
sx q[1];
rz(-1.5166538) q[1];
sx q[1];
rz(0.098310016) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5033103) q[0];
sx q[0];
rz(-1.1489963) q[0];
sx q[0];
rz(-0.038859239) q[0];
x q[1];
rz(3.1409522) q[2];
sx q[2];
rz(-1.5709236) q[2];
sx q[2];
rz(1.2544817) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6209539) q[1];
sx q[1];
rz(-2.6700182) q[1];
sx q[1];
rz(-3.0388799) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74048711) q[3];
sx q[3];
rz(-2.5038379) q[3];
sx q[3];
rz(-1.8000613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0564698) q[2];
sx q[2];
rz(-2.9462908) q[2];
sx q[2];
rz(-1.0657715) q[2];
rz(-0.39984518) q[3];
sx q[3];
rz(-2.6078434) q[3];
sx q[3];
rz(-1.2679509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9999076) q[0];
sx q[0];
rz(-0.081900224) q[0];
sx q[0];
rz(1.6837233) q[0];
rz(1.1204002) q[1];
sx q[1];
rz(-3.0034062) q[1];
sx q[1];
rz(-2.8021326) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2223245) q[0];
sx q[0];
rz(-1.647729) q[0];
sx q[0];
rz(-1.5841106) q[0];
rz(1.5846662) q[2];
sx q[2];
rz(-2.1254345) q[2];
sx q[2];
rz(-0.0019794606) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5684699) q[1];
sx q[1];
rz(-3.0559799) q[1];
sx q[1];
rz(1.5517615) q[1];
rz(-1.4711597) q[3];
sx q[3];
rz(-1.2335586) q[3];
sx q[3];
rz(1.3226313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7772943) q[2];
sx q[2];
rz(-0.0019145049) q[2];
sx q[2];
rz(0.36336362) q[2];
rz(-2.0511138) q[3];
sx q[3];
rz(-2.5616779) q[3];
sx q[3];
rz(-2.0232078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5217487) q[0];
sx q[0];
rz(-2.813297) q[0];
sx q[0];
rz(-2.2186665) q[0];
rz(-1.6587616) q[1];
sx q[1];
rz(-0.62983477) q[1];
sx q[1];
rz(3.1373851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.370458) q[0];
sx q[0];
rz(-2.1677289) q[0];
sx q[0];
rz(-2.202522) q[0];
rz(-pi) q[1];
rz(0.0059295456) q[2];
sx q[2];
rz(-1.1599132) q[2];
sx q[2];
rz(3.1290999) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9862822) q[1];
sx q[1];
rz(-2.661507) q[1];
sx q[1];
rz(-1.7106777) q[1];
rz(3.0757853) q[3];
sx q[3];
rz(-0.61344693) q[3];
sx q[3];
rz(2.891401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5620455) q[2];
sx q[2];
rz(-1.5374708) q[2];
sx q[2];
rz(1.9407678) q[2];
rz(-1.3935401) q[3];
sx q[3];
rz(-0.0037007185) q[3];
sx q[3];
rz(0.70011955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30018184) q[0];
sx q[0];
rz(-0.44898471) q[0];
sx q[0];
rz(1.9218943) q[0];
rz(1.3297184) q[1];
sx q[1];
rz(-1.1390353) q[1];
sx q[1];
rz(-0.16389287) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4224098) q[0];
sx q[0];
rz(-1.2162195) q[0];
sx q[0];
rz(0.18300458) q[0];
x q[1];
rz(0.40028769) q[2];
sx q[2];
rz(-0.62303715) q[2];
sx q[2];
rz(-1.3450587) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5245318) q[1];
sx q[1];
rz(-1.5111898) q[1];
sx q[1];
rz(-1.6888794) q[1];
rz(-pi) q[2];
rz(0.71299841) q[3];
sx q[3];
rz(-2.2621691) q[3];
sx q[3];
rz(2.612243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4280052) q[2];
sx q[2];
rz(-3.0516477) q[2];
sx q[2];
rz(-2.5186727) q[2];
rz(-0.04341393) q[3];
sx q[3];
rz(-0.89689887) q[3];
sx q[3];
rz(2.3626732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.6503705) q[0];
sx q[0];
rz(-3.1351884) q[0];
sx q[0];
rz(-2.6553335) q[0];
rz(-0.69475118) q[1];
sx q[1];
rz(-2.7943352) q[1];
sx q[1];
rz(2.6659226) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3846426) q[0];
sx q[0];
rz(-1.6044549) q[0];
sx q[0];
rz(-1.606694) q[0];
x q[1];
rz(1.6108405) q[2];
sx q[2];
rz(-2.4350428) q[2];
sx q[2];
rz(1.5020811) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6056532) q[1];
sx q[1];
rz(-1.5643143) q[1];
sx q[1];
rz(0.60757138) q[1];
x q[2];
rz(-0.31020152) q[3];
sx q[3];
rz(-1.7473012) q[3];
sx q[3];
rz(-1.4909397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.67705578) q[2];
sx q[2];
rz(-0.060364351) q[2];
sx q[2];
rz(-1.8439058) q[2];
rz(1.248598) q[3];
sx q[3];
rz(-2.5864351) q[3];
sx q[3];
rz(2.7738074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1260592) q[0];
sx q[0];
rz(-1.5801237) q[0];
sx q[0];
rz(1.7042241) q[0];
rz(2.2592648) q[1];
sx q[1];
rz(-0.066451646) q[1];
sx q[1];
rz(0.8393504) q[1];
rz(0.37226128) q[2];
sx q[2];
rz(-3.0426171) q[2];
sx q[2];
rz(-0.097886861) q[2];
rz(-0.23918693) q[3];
sx q[3];
rz(-1.5659955) q[3];
sx q[3];
rz(-1.5495054) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
