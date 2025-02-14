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
rz(1.7953402) q[0];
sx q[0];
rz(-1.6929545) q[0];
sx q[0];
rz(1.9899272) q[0];
rz(1.0898074) q[1];
sx q[1];
rz(-1.4767708) q[1];
sx q[1];
rz(-2.9906315) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38287791) q[0];
sx q[0];
rz(-2.8423514) q[0];
sx q[0];
rz(-1.3918357) q[0];
rz(-pi) q[1];
rz(1.483876) q[2];
sx q[2];
rz(-1.2582694) q[2];
sx q[2];
rz(-0.0258044) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.50569423) q[1];
sx q[1];
rz(-1.555996) q[1];
sx q[1];
rz(3.1379051) q[1];
x q[2];
rz(-1.5299705) q[3];
sx q[3];
rz(-1.4741952) q[3];
sx q[3];
rz(0.20303328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1208518) q[2];
sx q[2];
rz(-3.1181702) q[2];
sx q[2];
rz(-0.49129301) q[2];
rz(1.8190207) q[3];
sx q[3];
rz(-1.6072075) q[3];
sx q[3];
rz(1.8137431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604734) q[0];
sx q[0];
rz(-2.5871215) q[0];
sx q[0];
rz(2.4419899) q[0];
rz(-1.5910925) q[1];
sx q[1];
rz(-2.6467549) q[1];
sx q[1];
rz(2.9717305) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9014329) q[0];
sx q[0];
rz(-3.0548662) q[0];
sx q[0];
rz(-2.6804352) q[0];
rz(-0.40645941) q[2];
sx q[2];
rz(-0.32999295) q[2];
sx q[2];
rz(2.1292343) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.43812) q[1];
sx q[1];
rz(-1.2143163) q[1];
sx q[1];
rz(-2.0894024) q[1];
rz(-pi) q[2];
rz(2.027867) q[3];
sx q[3];
rz(-1.4997484) q[3];
sx q[3];
rz(0.20807264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.26213172) q[2];
sx q[2];
rz(-0.36236557) q[2];
sx q[2];
rz(1.7246838) q[2];
rz(1.0828241) q[3];
sx q[3];
rz(-2.2542451) q[3];
sx q[3];
rz(0.82720238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5812434) q[0];
sx q[0];
rz(-2.2214948) q[0];
sx q[0];
rz(-1.5592519) q[0];
rz(2.4371367) q[1];
sx q[1];
rz(-1.9943941) q[1];
sx q[1];
rz(-2.6047883) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0889595) q[0];
sx q[0];
rz(-0.68872243) q[0];
sx q[0];
rz(0.072865085) q[0];
x q[1];
rz(-0.090254668) q[2];
sx q[2];
rz(-1.5343416) q[2];
sx q[2];
rz(1.0807456) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4406261) q[1];
sx q[1];
rz(-1.4203209) q[1];
sx q[1];
rz(1.4626842) q[1];
x q[2];
rz(-2.267089) q[3];
sx q[3];
rz(-1.7603959) q[3];
sx q[3];
rz(-2.7087871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.392936) q[2];
sx q[2];
rz(-1.8042678) q[2];
sx q[2];
rz(-0.7788457) q[2];
rz(3.0574162) q[3];
sx q[3];
rz(-0.86936969) q[3];
sx q[3];
rz(-1.7387773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8970784) q[0];
sx q[0];
rz(-2.2205413) q[0];
sx q[0];
rz(2.670934) q[0];
rz(1.6828407) q[1];
sx q[1];
rz(-0.0048480821) q[1];
sx q[1];
rz(0.86476129) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4368468) q[0];
sx q[0];
rz(-0.517874) q[0];
sx q[0];
rz(0.11750607) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4385579) q[2];
sx q[2];
rz(-1.7900428) q[2];
sx q[2];
rz(-1.8869149) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4196849) q[1];
sx q[1];
rz(-3.0583755) q[1];
sx q[1];
rz(0.91099847) q[1];
x q[2];
rz(-1.8806609) q[3];
sx q[3];
rz(-1.052385) q[3];
sx q[3];
rz(-0.37004091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1137769) q[2];
sx q[2];
rz(-0.8646338) q[2];
sx q[2];
rz(1.9846385) q[2];
rz(2.786934) q[3];
sx q[3];
rz(-1.9054474) q[3];
sx q[3];
rz(-0.11053301) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.329634) q[0];
sx q[0];
rz(-2.2060153) q[0];
sx q[0];
rz(2.6079566) q[0];
rz(-1.2136906) q[1];
sx q[1];
rz(-3.1184989) q[1];
sx q[1];
rz(-1.1812706) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5438687) q[0];
sx q[0];
rz(-1.7337611) q[0];
sx q[0];
rz(-1.6173494) q[0];
x q[1];
rz(0.52631876) q[2];
sx q[2];
rz(-0.72229702) q[2];
sx q[2];
rz(2.1762951) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.54321982) q[1];
sx q[1];
rz(-0.66987619) q[1];
sx q[1];
rz(-1.3175439) q[1];
rz(-pi) q[2];
rz(0.13792928) q[3];
sx q[3];
rz(-1.0702881) q[3];
sx q[3];
rz(0.59964666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.10826762) q[2];
sx q[2];
rz(-2.6259618) q[2];
sx q[2];
rz(0.92833129) q[2];
rz(2.8632274) q[3];
sx q[3];
rz(-1.3185578) q[3];
sx q[3];
rz(0.068647169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5334897) q[0];
sx q[0];
rz(-2.6187496) q[0];
sx q[0];
rz(-0.19304092) q[0];
rz(2.4642384) q[1];
sx q[1];
rz(-3.1132128) q[1];
sx q[1];
rz(1.5100286) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0583147) q[0];
sx q[0];
rz(-1.4299576) q[0];
sx q[0];
rz(1.550564) q[0];
x q[1];
rz(-2.0329934) q[2];
sx q[2];
rz(-1.530556) q[2];
sx q[2];
rz(1.5239592) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2605156) q[1];
sx q[1];
rz(-2.279638) q[1];
sx q[1];
rz(-2.8167732) q[1];
x q[2];
rz(-0.72666177) q[3];
sx q[3];
rz(-1.9609986) q[3];
sx q[3];
rz(-1.778804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7789) q[2];
sx q[2];
rz(-0.71234667) q[2];
sx q[2];
rz(2.1544797) q[2];
rz(1.0846064) q[3];
sx q[3];
rz(-0.54607138) q[3];
sx q[3];
rz(0.62091056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5793594) q[0];
sx q[0];
rz(-2.1491304) q[0];
sx q[0];
rz(-1.5935422) q[0];
rz(0.69001895) q[1];
sx q[1];
rz(-3.1179805) q[1];
sx q[1];
rz(-1.4366368) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6520335) q[0];
sx q[0];
rz(-1.0512182) q[0];
sx q[0];
rz(1.621031) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68472548) q[2];
sx q[2];
rz(-1.5537571) q[2];
sx q[2];
rz(0.34340252) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1093398) q[1];
sx q[1];
rz(-2.5962862) q[1];
sx q[1];
rz(-1.4558353) q[1];
rz(-pi) q[2];
rz(1.8687181) q[3];
sx q[3];
rz(-2.0092589) q[3];
sx q[3];
rz(2.7496453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.37316608) q[2];
sx q[2];
rz(-1.9045786) q[2];
sx q[2];
rz(-0.58488673) q[2];
rz(-1.54555) q[3];
sx q[3];
rz(-1.7187748) q[3];
sx q[3];
rz(-0.91442951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.5717413) q[0];
sx q[0];
rz(-2.2241156) q[0];
sx q[0];
rz(1.5745987) q[0];
rz(0.51708108) q[1];
sx q[1];
rz(-2.2061429) q[1];
sx q[1];
rz(-1.2665952) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1388113) q[0];
sx q[0];
rz(-1.5673593) q[0];
sx q[0];
rz(3.1365738) q[0];
rz(2.3364303) q[2];
sx q[2];
rz(-1.1650208) q[2];
sx q[2];
rz(-2.9540887) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4609485) q[1];
sx q[1];
rz(-1.0310182) q[1];
sx q[1];
rz(2.7311027) q[1];
x q[2];
rz(0.1223337) q[3];
sx q[3];
rz(-1.4492826) q[3];
sx q[3];
rz(-2.7453855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7685585) q[2];
sx q[2];
rz(-2.3410102) q[2];
sx q[2];
rz(1.8689092) q[2];
rz(-2.7019165) q[3];
sx q[3];
rz(-1.9821854) q[3];
sx q[3];
rz(-3.0769949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30176297) q[0];
sx q[0];
rz(-0.01627144) q[0];
sx q[0];
rz(2.8453258) q[0];
rz(3.0773194) q[1];
sx q[1];
rz(-1.4646894) q[1];
sx q[1];
rz(-1.6305264) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26857527) q[0];
sx q[0];
rz(-2.7843256) q[0];
sx q[0];
rz(1.1910458) q[0];
rz(0.6728386) q[2];
sx q[2];
rz(-1.4110067) q[2];
sx q[2];
rz(0.21110134) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.293773) q[1];
sx q[1];
rz(-2.7917531) q[1];
sx q[1];
rz(-1.9560676) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1301163) q[3];
sx q[3];
rz(-1.6919129) q[3];
sx q[3];
rz(-0.83601421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54603798) q[2];
sx q[2];
rz(-1.7067355) q[2];
sx q[2];
rz(-1.1240553) q[2];
rz(1.3042287) q[3];
sx q[3];
rz(-0.29574695) q[3];
sx q[3];
rz(0.058163253) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8703363) q[0];
sx q[0];
rz(-0.79126343) q[0];
sx q[0];
rz(-1.1668209) q[0];
rz(-1.5406436) q[1];
sx q[1];
rz(-2.8438042) q[1];
sx q[1];
rz(1.3265532) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8542503) q[0];
sx q[0];
rz(-0.70579398) q[0];
sx q[0];
rz(-0.38278929) q[0];
x q[1];
rz(0.72188702) q[2];
sx q[2];
rz(-2.6195002) q[2];
sx q[2];
rz(2.4185857) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9611813) q[1];
sx q[1];
rz(-1.7827534) q[1];
sx q[1];
rz(2.4418264) q[1];
x q[2];
rz(-0.097685944) q[3];
sx q[3];
rz(-1.9703416) q[3];
sx q[3];
rz(-2.0942502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.88663179) q[2];
sx q[2];
rz(-2.9780264) q[2];
sx q[2];
rz(-1.4476042) q[2];
rz(-3.0149031) q[3];
sx q[3];
rz(-1.6619752) q[3];
sx q[3];
rz(1.1160342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.806725) q[0];
sx q[0];
rz(-1.3074449) q[0];
sx q[0];
rz(1.773651) q[0];
rz(-1.5927636) q[1];
sx q[1];
rz(-0.82850414) q[1];
sx q[1];
rz(-2.9864476) q[1];
rz(-2.7663076) q[2];
sx q[2];
rz(-0.47561687) q[2];
sx q[2];
rz(0.33676666) q[2];
rz(-2.9709076) q[3];
sx q[3];
rz(-0.99320929) q[3];
sx q[3];
rz(0.39218642) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
