OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(-3.0769899) q[0];
sx q[0];
rz(-0.021615418) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(1.8145476) q[1];
sx q[1];
rz(10.756455) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0639122) q[0];
sx q[0];
rz(-2.4644116) q[0];
sx q[0];
rz(2.1190686) q[0];
rz(-pi) q[1];
rz(0.10279074) q[2];
sx q[2];
rz(-1.79709) q[2];
sx q[2];
rz(-3.0068827) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.057030765) q[1];
sx q[1];
rz(-1.7006526) q[1];
sx q[1];
rz(-0.76743857) q[1];
rz(-pi) q[2];
rz(0.33820037) q[3];
sx q[3];
rz(-2.1669398) q[3];
sx q[3];
rz(0.81830922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8109479) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(-0.60418207) q[2];
rz(-1.0243246) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(2.0143051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3246831) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(-2.0781793) q[0];
rz(2.2564607) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(3.1399472) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9379399) q[0];
sx q[0];
rz(-1.5814039) q[0];
sx q[0];
rz(-1.525735) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64237853) q[2];
sx q[2];
rz(-1.3628236) q[2];
sx q[2];
rz(0.013052879) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8811223) q[1];
sx q[1];
rz(-2.3350041) q[1];
sx q[1];
rz(1.9116198) q[1];
x q[2];
rz(0.1699154) q[3];
sx q[3];
rz(-1.8222457) q[3];
sx q[3];
rz(2.2268695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1559747) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(0.70811159) q[2];
rz(-2.0478785) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(0.99350199) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.920632) q[0];
sx q[0];
rz(-1.4039803) q[0];
sx q[0];
rz(0.8272585) q[0];
rz(-0.0050841252) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(-1.089383) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2978965) q[0];
sx q[0];
rz(-1.1481522) q[0];
sx q[0];
rz(-2.1556913) q[0];
rz(-pi) q[1];
rz(-0.59805869) q[2];
sx q[2];
rz(-1.2130249) q[2];
sx q[2];
rz(1.9269112) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1091724) q[1];
sx q[1];
rz(-1.4854327) q[1];
sx q[1];
rz(-1.1847772) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.009922) q[3];
sx q[3];
rz(-0.86285931) q[3];
sx q[3];
rz(2.5885955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0638782) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(2.2568978) q[2];
rz(1.9583154) q[3];
sx q[3];
rz(-1.3180472) q[3];
sx q[3];
rz(-1.6515091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5723715) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(-2.4147721) q[0];
rz(-0.80365333) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(2.7817536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026304631) q[0];
sx q[0];
rz(-1.7548772) q[0];
sx q[0];
rz(-2.227965) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6511351) q[2];
sx q[2];
rz(-1.6614117) q[2];
sx q[2];
rz(-0.091094253) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.53107809) q[1];
sx q[1];
rz(-0.58137608) q[1];
sx q[1];
rz(-2.4060712) q[1];
x q[2];
rz(1.2218277) q[3];
sx q[3];
rz(-1.2811321) q[3];
sx q[3];
rz(1.3834013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.56746733) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(2.3763669) q[2];
rz(2.3848173) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(-0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9969479) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(1.7255406) q[0];
rz(-2.7744746) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(-1.0353154) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9177592) q[0];
sx q[0];
rz(-2.1609554) q[0];
sx q[0];
rz(0.53442861) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5315227) q[2];
sx q[2];
rz(-2.0619259) q[2];
sx q[2];
rz(2.479535) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.94381911) q[1];
sx q[1];
rz(-0.10995956) q[1];
sx q[1];
rz(-2.6206559) q[1];
rz(-pi) q[2];
rz(-2.4953793) q[3];
sx q[3];
rz(-0.34020243) q[3];
sx q[3];
rz(0.37341213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3570024) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(-0.76888293) q[2];
rz(0.33603493) q[3];
sx q[3];
rz(-2.3612645) q[3];
sx q[3];
rz(1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96200213) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(1.996421) q[0];
rz(-1.1046474) q[1];
sx q[1];
rz(-1.2303753) q[1];
sx q[1];
rz(-0.2072269) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4743054) q[0];
sx q[0];
rz(-1.4315839) q[0];
sx q[0];
rz(2.5179203) q[0];
x q[1];
rz(-0.63164288) q[2];
sx q[2];
rz(-2.4372059) q[2];
sx q[2];
rz(-0.44142516) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6991899) q[1];
sx q[1];
rz(-1.8132205) q[1];
sx q[1];
rz(-1.8276023) q[1];
rz(1.7436142) q[3];
sx q[3];
rz(-2.380905) q[3];
sx q[3];
rz(-2.2067604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8032288) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(0.72171372) q[2];
rz(-1.2747814) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(2.8815564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.24467829) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(2.4095643) q[0];
rz(-0.0094982068) q[1];
sx q[1];
rz(-0.54832012) q[1];
sx q[1];
rz(2.8498555) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1222629) q[0];
sx q[0];
rz(-2.7195192) q[0];
sx q[0];
rz(-3.0164369) q[0];
rz(1.8998434) q[2];
sx q[2];
rz(-2.8223158) q[2];
sx q[2];
rz(0.74908756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2597255) q[1];
sx q[1];
rz(-1.423466) q[1];
sx q[1];
rz(-2.8357361) q[1];
rz(0.52280207) q[3];
sx q[3];
rz(-1.3005946) q[3];
sx q[3];
rz(3.0772046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26414028) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(-2.3042802) q[2];
rz(-1.1710179) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(3.0795735) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0432805) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(-0.6638546) q[0];
rz(3.035416) q[1];
sx q[1];
rz(-2.5352434) q[1];
sx q[1];
rz(0.95867872) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6129235) q[0];
sx q[0];
rz(-1.0597502) q[0];
sx q[0];
rz(0.59707609) q[0];
rz(0.27997048) q[2];
sx q[2];
rz(-1.0500056) q[2];
sx q[2];
rz(2.1998646) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.16756032) q[1];
sx q[1];
rz(-1.6775727) q[1];
sx q[1];
rz(-2.0978931) q[1];
x q[2];
rz(-1.6931375) q[3];
sx q[3];
rz(-0.99007505) q[3];
sx q[3];
rz(1.3337222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.133193) q[2];
sx q[2];
rz(-1.2352751) q[2];
sx q[2];
rz(-0.91910249) q[2];
rz(-1.3778) q[3];
sx q[3];
rz(-1.931124) q[3];
sx q[3];
rz(1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.33525) q[0];
sx q[0];
rz(-0.85334539) q[0];
sx q[0];
rz(-2.7046955) q[0];
rz(0.70029744) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(2.6760496) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7882999) q[0];
sx q[0];
rz(-0.71209891) q[0];
sx q[0];
rz(-3.0112991) q[0];
rz(-pi) q[1];
rz(3.0131857) q[2];
sx q[2];
rz(-1.1616716) q[2];
sx q[2];
rz(-0.94305925) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1381582) q[1];
sx q[1];
rz(-2.5596566) q[1];
sx q[1];
rz(2.483063) q[1];
x q[2];
rz(-2.8903928) q[3];
sx q[3];
rz(-1.9198872) q[3];
sx q[3];
rz(-0.72124764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(-1.4979866) q[2];
rz(2.9368029) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(-2.8695316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64514226) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(-2.0196594) q[0];
rz(-2.3902068) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(0.5823935) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5781007) q[0];
sx q[0];
rz(-2.2390215) q[0];
sx q[0];
rz(0.25989728) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15162823) q[2];
sx q[2];
rz(-2.0564338) q[2];
sx q[2];
rz(-3.1079353) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45422428) q[1];
sx q[1];
rz(-2.9240989) q[1];
sx q[1];
rz(-0.99517676) q[1];
rz(-0.12606975) q[3];
sx q[3];
rz(-0.21890103) q[3];
sx q[3];
rz(-2.253988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0104684) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(-2.3790512) q[2];
rz(1.4108346) q[3];
sx q[3];
rz(-2.2236731) q[3];
sx q[3];
rz(-1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0653771) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(-1.8021884) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(0.28152485) q[2];
sx q[2];
rz(-0.83419656) q[2];
sx q[2];
rz(0.4783334) q[2];
rz(-2.6639832) q[3];
sx q[3];
rz(-1.1160679) q[3];
sx q[3];
rz(3.06649) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
