OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.66981411) q[0];
sx q[0];
rz(1.9999003) q[0];
sx q[0];
rz(6.5394149) q[0];
rz(0.37707075) q[1];
sx q[1];
rz(-1.7989676) q[1];
sx q[1];
rz(-2.0816198) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.357517) q[0];
sx q[0];
rz(-1.6889484) q[0];
sx q[0];
rz(-1.4489001) q[0];
rz(-pi) q[1];
rz(0.04214123) q[2];
sx q[2];
rz(-2.4486447) q[2];
sx q[2];
rz(2.9646089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4958475) q[1];
sx q[1];
rz(-1.8105257) q[1];
sx q[1];
rz(1.354753) q[1];
rz(-2.5422402) q[3];
sx q[3];
rz(-2.9671892) q[3];
sx q[3];
rz(2.6389183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18067351) q[2];
sx q[2];
rz(-2.2283165) q[2];
sx q[2];
rz(1.712435) q[2];
rz(-0.6116496) q[3];
sx q[3];
rz(-1.4432171) q[3];
sx q[3];
rz(3.070224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3058158) q[0];
sx q[0];
rz(-2.340305) q[0];
sx q[0];
rz(2.3961156) q[0];
rz(1.9396793) q[1];
sx q[1];
rz(-1.3246091) q[1];
sx q[1];
rz(-1.036693) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1069477) q[0];
sx q[0];
rz(-1.6634403) q[0];
sx q[0];
rz(-1.9361467) q[0];
x q[1];
rz(3.1142508) q[2];
sx q[2];
rz(-0.76702709) q[2];
sx q[2];
rz(2.1576144) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9696641) q[1];
sx q[1];
rz(-1.2873889) q[1];
sx q[1];
rz(-0.71782459) q[1];
x q[2];
rz(-0.58425957) q[3];
sx q[3];
rz(-1.7528894) q[3];
sx q[3];
rz(0.77891536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1408954) q[2];
sx q[2];
rz(-1.8706198) q[2];
sx q[2];
rz(2.1689283) q[2];
rz(-2.2843212) q[3];
sx q[3];
rz(-1.0664777) q[3];
sx q[3];
rz(1.8734141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5001517) q[0];
sx q[0];
rz(-0.86857906) q[0];
sx q[0];
rz(0.80297536) q[0];
rz(0.650644) q[1];
sx q[1];
rz(-2.3425808) q[1];
sx q[1];
rz(-2.9237936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2265711) q[0];
sx q[0];
rz(-2.3816097) q[0];
sx q[0];
rz(2.1334011) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62816633) q[2];
sx q[2];
rz(-1.7587314) q[2];
sx q[2];
rz(-2.4357901) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7347225) q[1];
sx q[1];
rz(-1.0823421) q[1];
sx q[1];
rz(0.10481701) q[1];
rz(2.5434615) q[3];
sx q[3];
rz(-0.39970988) q[3];
sx q[3];
rz(-2.3593115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.95152068) q[2];
sx q[2];
rz(-1.051544) q[2];
sx q[2];
rz(-1.5104843) q[2];
rz(-1.914628) q[3];
sx q[3];
rz(-2.0881784) q[3];
sx q[3];
rz(2.1372883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.889582) q[0];
sx q[0];
rz(-2.4376106) q[0];
sx q[0];
rz(-0.66396436) q[0];
rz(-0.12531677) q[1];
sx q[1];
rz(-1.434451) q[1];
sx q[1];
rz(-1.0391611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9665899) q[0];
sx q[0];
rz(-3.1313854) q[0];
sx q[0];
rz(-0.27423476) q[0];
rz(2.152056) q[2];
sx q[2];
rz(-0.46316812) q[2];
sx q[2];
rz(0.65639596) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3287828) q[1];
sx q[1];
rz(-0.53063466) q[1];
sx q[1];
rz(-1.1483436) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6986946) q[3];
sx q[3];
rz(-1.6567461) q[3];
sx q[3];
rz(-1.3493054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8404954) q[2];
sx q[2];
rz(-1.4157462) q[2];
sx q[2];
rz(-0.077979716) q[2];
rz(1.1178499) q[3];
sx q[3];
rz(-1.0406787) q[3];
sx q[3];
rz(-1.0361205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2074821) q[0];
sx q[0];
rz(-0.20620646) q[0];
sx q[0];
rz(2.7284486) q[0];
rz(1.9059937) q[1];
sx q[1];
rz(-1.8159591) q[1];
sx q[1];
rz(1.3135757) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51172709) q[0];
sx q[0];
rz(-2.4851375) q[0];
sx q[0];
rz(1.8740428) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1059472) q[2];
sx q[2];
rz(-1.2463971) q[2];
sx q[2];
rz(1.731763) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8660248) q[1];
sx q[1];
rz(-1.1236262) q[1];
sx q[1];
rz(-0.42110301) q[1];
rz(-pi) q[2];
rz(2.3543139) q[3];
sx q[3];
rz(-1.9036999) q[3];
sx q[3];
rz(1.5008139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1818587) q[2];
sx q[2];
rz(-1.2343957) q[2];
sx q[2];
rz(2.6971297) q[2];
rz(1.2005165) q[3];
sx q[3];
rz(-1.8903172) q[3];
sx q[3];
rz(-0.1058696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7926517) q[0];
sx q[0];
rz(-2.7466725) q[0];
sx q[0];
rz(-3.0643903) q[0];
rz(0.96241799) q[1];
sx q[1];
rz(-1.9541157) q[1];
sx q[1];
rz(-0.033800689) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156436) q[0];
sx q[0];
rz(-1.9140049) q[0];
sx q[0];
rz(2.6604466) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5438429) q[2];
sx q[2];
rz(-1.0232864) q[2];
sx q[2];
rz(-1.6545878) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5742366) q[1];
sx q[1];
rz(-2.1891382) q[1];
sx q[1];
rz(-0.38783698) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.364879) q[3];
sx q[3];
rz(-1.6459822) q[3];
sx q[3];
rz(3.1326339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9683711) q[2];
sx q[2];
rz(-1.5366448) q[2];
sx q[2];
rz(2.8823493) q[2];
rz(0.94414583) q[3];
sx q[3];
rz(-0.21374948) q[3];
sx q[3];
rz(-1.3313782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060870085) q[0];
sx q[0];
rz(-0.20651564) q[0];
sx q[0];
rz(2.5282705) q[0];
rz(-2.2382286) q[1];
sx q[1];
rz(-0.84066835) q[1];
sx q[1];
rz(0.14796743) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6981993) q[0];
sx q[0];
rz(-0.86981378) q[0];
sx q[0];
rz(-2.2156391) q[0];
x q[1];
rz(1.5388707) q[2];
sx q[2];
rz(-2.0096001) q[2];
sx q[2];
rz(1.9500458) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.24616775) q[1];
sx q[1];
rz(-2.2895799) q[1];
sx q[1];
rz(-2.5661929) q[1];
rz(-pi) q[2];
rz(2.1359753) q[3];
sx q[3];
rz(-1.5907173) q[3];
sx q[3];
rz(-1.3764149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2030187) q[2];
sx q[2];
rz(-1.4116986) q[2];
sx q[2];
rz(-0.99986783) q[2];
rz(1.8262919) q[3];
sx q[3];
rz(-0.2838997) q[3];
sx q[3];
rz(2.4711173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1891747) q[0];
sx q[0];
rz(-0.84680951) q[0];
sx q[0];
rz(2.1536105) q[0];
rz(1.2692163) q[1];
sx q[1];
rz(-0.80943426) q[1];
sx q[1];
rz(2.9775137) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39597797) q[0];
sx q[0];
rz(-1.4854637) q[0];
sx q[0];
rz(3.0239168) q[0];
rz(-pi) q[1];
rz(2.9226926) q[2];
sx q[2];
rz(-2.0563807) q[2];
sx q[2];
rz(1.4771493) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.33262256) q[1];
sx q[1];
rz(-1.5676985) q[1];
sx q[1];
rz(1.078152) q[1];
rz(-pi) q[2];
rz(2.2211391) q[3];
sx q[3];
rz(-1.4152539) q[3];
sx q[3];
rz(0.33442861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24171955) q[2];
sx q[2];
rz(-2.6417929) q[2];
sx q[2];
rz(-1.0618173) q[2];
rz(-3.0070983) q[3];
sx q[3];
rz(-0.89956346) q[3];
sx q[3];
rz(-3.0883446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4835085) q[0];
sx q[0];
rz(-0.71165076) q[0];
sx q[0];
rz(-0.2051951) q[0];
rz(-0.23458734) q[1];
sx q[1];
rz(-1.4261475) q[1];
sx q[1];
rz(-0.1964143) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097986624) q[0];
sx q[0];
rz(-2.4825486) q[0];
sx q[0];
rz(-2.090467) q[0];
rz(-pi) q[1];
rz(1.4142198) q[2];
sx q[2];
rz(-1.051826) q[2];
sx q[2];
rz(-1.165498) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3383693) q[1];
sx q[1];
rz(-1.3937989) q[1];
sx q[1];
rz(-1.2425674) q[1];
rz(-pi) q[2];
rz(-0.78077448) q[3];
sx q[3];
rz(-1.1765683) q[3];
sx q[3];
rz(0.0782644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3966169) q[2];
sx q[2];
rz(-1.1810415) q[2];
sx q[2];
rz(1.1783925) q[2];
rz(0.34934238) q[3];
sx q[3];
rz(-1.1290461) q[3];
sx q[3];
rz(1.0497302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.0013393764) q[0];
sx q[0];
rz(-1.1118735) q[0];
sx q[0];
rz(-2.4358791) q[0];
rz(-0.69752518) q[1];
sx q[1];
rz(-1.3730647) q[1];
sx q[1];
rz(-2.2360905) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4587935) q[0];
sx q[0];
rz(-2.0044649) q[0];
sx q[0];
rz(-1.5894029) q[0];
x q[1];
rz(1.2442676) q[2];
sx q[2];
rz(-1.2249607) q[2];
sx q[2];
rz(-2.2188733) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6433555) q[1];
sx q[1];
rz(-1.61803) q[1];
sx q[1];
rz(2.8308353) q[1];
x q[2];
rz(1.2071183) q[3];
sx q[3];
rz(-1.8810026) q[3];
sx q[3];
rz(0.23387533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8991578) q[2];
sx q[2];
rz(-2.0826715) q[2];
sx q[2];
rz(-2.4230797) q[2];
rz(2.4588623) q[3];
sx q[3];
rz(-1.9971137) q[3];
sx q[3];
rz(-3.0375286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.318442) q[0];
sx q[0];
rz(-2.5704076) q[0];
sx q[0];
rz(3.0388863) q[0];
rz(1.6406583) q[1];
sx q[1];
rz(-1.8713015) q[1];
sx q[1];
rz(-0.4458977) q[1];
rz(1.4529543) q[2];
sx q[2];
rz(-1.4284882) q[2];
sx q[2];
rz(1.535648) q[2];
rz(-2.2352696) q[3];
sx q[3];
rz(-0.88574468) q[3];
sx q[3];
rz(1.1024324) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
