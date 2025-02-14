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
rz(2.6858202) q[0];
sx q[0];
rz(-2.0210285) q[0];
sx q[0];
rz(1.7280581) q[0];
rz(-2.7780374) q[1];
sx q[1];
rz(-1.8283365) q[1];
sx q[1];
rz(0.29120905) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7257627) q[0];
sx q[0];
rz(-2.8208469) q[0];
sx q[0];
rz(-1.6345535) q[0];
x q[1];
rz(-1.2313263) q[2];
sx q[2];
rz(-1.4541153) q[2];
sx q[2];
rz(1.1543903) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26303459) q[1];
sx q[1];
rz(-2.546995) q[1];
sx q[1];
rz(0.27185284) q[1];
rz(-pi) q[2];
rz(-2.5194571) q[3];
sx q[3];
rz(-0.91968007) q[3];
sx q[3];
rz(1.5762005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1353961) q[2];
sx q[2];
rz(-1.5781032) q[2];
sx q[2];
rz(-0.023999365) q[2];
rz(0.35605797) q[3];
sx q[3];
rz(-0.67432109) q[3];
sx q[3];
rz(0.02478987) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0559219) q[0];
sx q[0];
rz(-2.4503158) q[0];
sx q[0];
rz(0.63496494) q[0];
rz(-0.7629281) q[1];
sx q[1];
rz(-1.678391) q[1];
sx q[1];
rz(-0.91711226) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35058366) q[0];
sx q[0];
rz(-1.6054594) q[0];
sx q[0];
rz(-3.123218) q[0];
rz(2.2516052) q[2];
sx q[2];
rz(-1.1900848) q[2];
sx q[2];
rz(3.0842881) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8626111) q[1];
sx q[1];
rz(-1.647842) q[1];
sx q[1];
rz(2.0120697) q[1];
x q[2];
rz(-0.95083745) q[3];
sx q[3];
rz(-1.3966433) q[3];
sx q[3];
rz(0.0096498409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1941173) q[2];
sx q[2];
rz(-0.52299356) q[2];
sx q[2];
rz(0.73327649) q[2];
rz(-2.0745847) q[3];
sx q[3];
rz(-1.7206444) q[3];
sx q[3];
rz(-0.1799306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5718403) q[0];
sx q[0];
rz(-1.1886007) q[0];
sx q[0];
rz(-2.6133614) q[0];
rz(0.86604467) q[1];
sx q[1];
rz(-1.8128315) q[1];
sx q[1];
rz(-2.6285062) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24241867) q[0];
sx q[0];
rz(-2.0089825) q[0];
sx q[0];
rz(-0.04695453) q[0];
x q[1];
rz(0.92136356) q[2];
sx q[2];
rz(-1.3201534) q[2];
sx q[2];
rz(-1.7684557) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5410903) q[1];
sx q[1];
rz(-1.2230495) q[1];
sx q[1];
rz(0.13558023) q[1];
rz(-pi) q[2];
rz(-3.0894075) q[3];
sx q[3];
rz(-2.6413249) q[3];
sx q[3];
rz(-2.1310998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1055866) q[2];
sx q[2];
rz(-1.0893818) q[2];
sx q[2];
rz(3.1028683) q[2];
rz(0.45768467) q[3];
sx q[3];
rz(-0.80515146) q[3];
sx q[3];
rz(-1.800644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4337346) q[0];
sx q[0];
rz(-2.7271294) q[0];
sx q[0];
rz(1.1210972) q[0];
rz(2.3838249) q[1];
sx q[1];
rz(-1.517375) q[1];
sx q[1];
rz(-0.68414348) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093903784) q[0];
sx q[0];
rz(-2.2950116) q[0];
sx q[0];
rz(2.8709477) q[0];
rz(-2.5252417) q[2];
sx q[2];
rz(-3.1292097) q[2];
sx q[2];
rz(-1.826926) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.63057478) q[1];
sx q[1];
rz(-1.2666432) q[1];
sx q[1];
rz(-1.263878) q[1];
rz(-pi) q[2];
rz(-2.6287967) q[3];
sx q[3];
rz(-1.0362451) q[3];
sx q[3];
rz(-3.0844546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6571558) q[2];
sx q[2];
rz(-1.6146722) q[2];
sx q[2];
rz(0.26051513) q[2];
rz(-2.303463) q[3];
sx q[3];
rz(-1.9939491) q[3];
sx q[3];
rz(-0.014613541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92358661) q[0];
sx q[0];
rz(-1.7638693) q[0];
sx q[0];
rz(-2.6309784) q[0];
rz(1.8922197) q[1];
sx q[1];
rz(-1.1621954) q[1];
sx q[1];
rz(-1.6872663) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73390612) q[0];
sx q[0];
rz(-2.9260194) q[0];
sx q[0];
rz(0.13224299) q[0];
rz(-pi) q[1];
rz(-3.0018423) q[2];
sx q[2];
rz(-1.5232161) q[2];
sx q[2];
rz(-0.51932491) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1087423) q[1];
sx q[1];
rz(-1.3128442) q[1];
sx q[1];
rz(-0.56405073) q[1];
x q[2];
rz(0.14438676) q[3];
sx q[3];
rz(-1.4098415) q[3];
sx q[3];
rz(-2.150589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.55383596) q[2];
sx q[2];
rz(-0.38148701) q[2];
sx q[2];
rz(2.7510263) q[2];
rz(-0.25911123) q[3];
sx q[3];
rz(-1.6389537) q[3];
sx q[3];
rz(-0.34671569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7873586) q[0];
sx q[0];
rz(-2.7095095) q[0];
sx q[0];
rz(-0.64467347) q[0];
rz(-1.3020172) q[1];
sx q[1];
rz(-1.6141067) q[1];
sx q[1];
rz(-2.7929746) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7946588) q[0];
sx q[0];
rz(-1.1125135) q[0];
sx q[0];
rz(0.51308227) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.604282) q[2];
sx q[2];
rz(-1.9331353) q[2];
sx q[2];
rz(0.57149959) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4069566) q[1];
sx q[1];
rz(-0.52328456) q[1];
sx q[1];
rz(0.95013817) q[1];
rz(-1.5362605) q[3];
sx q[3];
rz(-1.4538962) q[3];
sx q[3];
rz(-1.4394606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0384486) q[2];
sx q[2];
rz(-0.95547533) q[2];
sx q[2];
rz(-0.038334282) q[2];
rz(-0.96758715) q[3];
sx q[3];
rz(-1.4097694) q[3];
sx q[3];
rz(2.7819395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86647812) q[0];
sx q[0];
rz(-2.622128) q[0];
sx q[0];
rz(3.0679829) q[0];
rz(1.7346802) q[1];
sx q[1];
rz(-0.55714566) q[1];
sx q[1];
rz(0.73307347) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063696472) q[0];
sx q[0];
rz(-1.0390345) q[0];
sx q[0];
rz(-0.97490208) q[0];
rz(2.2233811) q[2];
sx q[2];
rz(-0.84647734) q[2];
sx q[2];
rz(2.0744963) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9736874) q[1];
sx q[1];
rz(-1.2406772) q[1];
sx q[1];
rz(-2.7658505) q[1];
x q[2];
rz(-0.2352318) q[3];
sx q[3];
rz(-1.8149655) q[3];
sx q[3];
rz(-1.5489123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0343895) q[2];
sx q[2];
rz(-1.723039) q[2];
sx q[2];
rz(3.1305195) q[2];
rz(2.8360046) q[3];
sx q[3];
rz(-1.1253076) q[3];
sx q[3];
rz(-1.2529469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3333862) q[0];
sx q[0];
rz(-0.60788637) q[0];
sx q[0];
rz(-0.73390865) q[0];
rz(-0.58087307) q[1];
sx q[1];
rz(-1.0092694) q[1];
sx q[1];
rz(3.0427921) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0990813) q[0];
sx q[0];
rz(-2.3020199) q[0];
sx q[0];
rz(0.73584105) q[0];
x q[1];
rz(-0.96510624) q[2];
sx q[2];
rz(-0.98585502) q[2];
sx q[2];
rz(-3.0966058) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.172625) q[1];
sx q[1];
rz(-1.8207112) q[1];
sx q[1];
rz(2.6843798) q[1];
rz(-pi) q[2];
rz(0.55136724) q[3];
sx q[3];
rz(-0.15227642) q[3];
sx q[3];
rz(0.65413977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8431479) q[2];
sx q[2];
rz(-1.2848022) q[2];
sx q[2];
rz(1.0996381) q[2];
rz(0.034865033) q[3];
sx q[3];
rz(-2.3979082) q[3];
sx q[3];
rz(0.59665027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027503969) q[0];
sx q[0];
rz(-1.9845668) q[0];
sx q[0];
rz(-1.9644568) q[0];
rz(-1.4741395) q[1];
sx q[1];
rz(-2.2087704) q[1];
sx q[1];
rz(1.4449545) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.225017) q[0];
sx q[0];
rz(-1.9635884) q[0];
sx q[0];
rz(1.265822) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99924008) q[2];
sx q[2];
rz(-0.88724595) q[2];
sx q[2];
rz(0.37887606) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8813546) q[1];
sx q[1];
rz(-1.5020276) q[1];
sx q[1];
rz(-2.0860904) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8650896) q[3];
sx q[3];
rz(-0.91405896) q[3];
sx q[3];
rz(1.2402676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16704796) q[2];
sx q[2];
rz(-1.5844774) q[2];
sx q[2];
rz(0.17189279) q[2];
rz(-0.78957549) q[3];
sx q[3];
rz(-1.1372477) q[3];
sx q[3];
rz(1.3353698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0563141) q[0];
sx q[0];
rz(-2.9264937) q[0];
sx q[0];
rz(-2.6000182) q[0];
rz(-1.1594634) q[1];
sx q[1];
rz(-1.3075202) q[1];
sx q[1];
rz(-0.83311876) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3726927) q[0];
sx q[0];
rz(-1.8463377) q[0];
sx q[0];
rz(-2.5252456) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6037349) q[2];
sx q[2];
rz(-0.61213697) q[2];
sx q[2];
rz(-2.9878841) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45634746) q[1];
sx q[1];
rz(-0.73328555) q[1];
sx q[1];
rz(1.9128146) q[1];
x q[2];
rz(-1.5686137) q[3];
sx q[3];
rz(-2.4139521) q[3];
sx q[3];
rz(-1.997987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1200166) q[2];
sx q[2];
rz(-2.3638201) q[2];
sx q[2];
rz(-0.94338256) q[2];
rz(-1.0776445) q[3];
sx q[3];
rz(-1.5871983) q[3];
sx q[3];
rz(2.8130892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42644603) q[0];
sx q[0];
rz(-1.0363415) q[0];
sx q[0];
rz(2.8847726) q[0];
rz(-0.23964755) q[1];
sx q[1];
rz(-0.9050723) q[1];
sx q[1];
rz(-1.1368652) q[1];
rz(0.69375615) q[2];
sx q[2];
rz(-2.0907331) q[2];
sx q[2];
rz(1.5246684) q[2];
rz(1.3171468) q[3];
sx q[3];
rz(-1.8719049) q[3];
sx q[3];
rz(-0.20897839) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
