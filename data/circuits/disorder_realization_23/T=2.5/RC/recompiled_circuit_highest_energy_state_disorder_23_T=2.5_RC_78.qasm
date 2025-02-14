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
rz(-1.3462525) q[0];
sx q[0];
rz(-1.4486382) q[0];
sx q[0];
rz(1.1516655) q[0];
rz(1.0898074) q[1];
sx q[1];
rz(-1.4767708) q[1];
sx q[1];
rz(-2.9906315) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7587147) q[0];
sx q[0];
rz(-0.29924127) q[0];
sx q[0];
rz(1.749757) q[0];
x q[1];
rz(-0.26246983) q[2];
sx q[2];
rz(-2.8175857) q[2];
sx q[2];
rz(-2.8912247) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0651567) q[1];
sx q[1];
rz(-1.5744835) q[1];
sx q[1];
rz(1.5855968) q[1];
x q[2];
rz(-0.096681194) q[3];
sx q[3];
rz(-1.6114317) q[3];
sx q[3];
rz(1.3717029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0207409) q[2];
sx q[2];
rz(-0.023422478) q[2];
sx q[2];
rz(-0.49129301) q[2];
rz(-1.8190207) q[3];
sx q[3];
rz(-1.6072075) q[3];
sx q[3];
rz(1.3278495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604734) q[0];
sx q[0];
rz(-2.5871215) q[0];
sx q[0];
rz(-2.4419899) q[0];
rz(1.5505002) q[1];
sx q[1];
rz(-0.49483776) q[1];
sx q[1];
rz(0.16986212) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79029523) q[0];
sx q[0];
rz(-1.6093495) q[0];
sx q[0];
rz(0.077705381) q[0];
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
rz(-0.58386034) q[1];
sx q[1];
rz(-0.61994821) q[1];
sx q[1];
rz(2.2151193) q[1];
rz(-pi) q[2];
rz(1.1137257) q[3];
sx q[3];
rz(-1.6418442) q[3];
sx q[3];
rz(-2.93352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8794609) q[2];
sx q[2];
rz(-2.7792271) q[2];
sx q[2];
rz(1.4169089) q[2];
rz(1.0828241) q[3];
sx q[3];
rz(-2.2542451) q[3];
sx q[3];
rz(0.82720238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5603492) q[0];
sx q[0];
rz(-2.2214948) q[0];
sx q[0];
rz(1.5823407) q[0];
rz(2.4371367) q[1];
sx q[1];
rz(-1.1471986) q[1];
sx q[1];
rz(2.6047883) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0999609) q[0];
sx q[0];
rz(-2.2573364) q[0];
sx q[0];
rz(-1.6306535) q[0];
x q[1];
rz(0.38449826) q[2];
sx q[2];
rz(-0.097320312) q[2];
sx q[2];
rz(-2.26869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0672863) q[1];
sx q[1];
rz(-0.18504772) q[1];
sx q[1];
rz(0.61850278) q[1];
x q[2];
rz(-1.8615253) q[3];
sx q[3];
rz(-0.71746263) q[3];
sx q[3];
rz(0.91625253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.74865666) q[2];
sx q[2];
rz(-1.8042678) q[2];
sx q[2];
rz(-0.7788457) q[2];
rz(-3.0574162) q[3];
sx q[3];
rz(-0.86936969) q[3];
sx q[3];
rz(1.7387773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.8970784) q[0];
sx q[0];
rz(-0.92105138) q[0];
sx q[0];
rz(0.47065863) q[0];
rz(-1.6828407) q[1];
sx q[1];
rz(-3.1367446) q[1];
sx q[1];
rz(-2.2768314) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1098561) q[0];
sx q[0];
rz(-1.512728) q[0];
sx q[0];
rz(-0.51490291) q[0];
rz(-pi) q[1];
rz(0.22111544) q[2];
sx q[2];
rz(-1.4417414) q[2];
sx q[2];
rz(-0.2871967) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7219077) q[1];
sx q[1];
rz(-0.083217155) q[1];
sx q[1];
rz(-2.2305942) q[1];
rz(-pi) q[2];
rz(-1.8806609) q[3];
sx q[3];
rz(-2.0892077) q[3];
sx q[3];
rz(-2.7715517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1137769) q[2];
sx q[2];
rz(-2.2769589) q[2];
sx q[2];
rz(-1.9846385) q[2];
rz(-2.786934) q[3];
sx q[3];
rz(-1.2361453) q[3];
sx q[3];
rz(3.0310596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.329634) q[0];
sx q[0];
rz(-0.93557731) q[0];
sx q[0];
rz(-0.53363609) q[0];
rz(1.927902) q[1];
sx q[1];
rz(-0.023093725) q[1];
sx q[1];
rz(1.1812706) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.597724) q[0];
sx q[0];
rz(-1.7337611) q[0];
sx q[0];
rz(-1.5242432) q[0];
rz(-1.9875162) q[2];
sx q[2];
rz(-0.96229711) q[2];
sx q[2];
rz(1.5173943) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.279484) q[1];
sx q[1];
rz(-0.92593926) q[1];
sx q[1];
rz(-0.1959066) q[1];
rz(-pi) q[2];
x q[2];
rz(1.817068) q[3];
sx q[3];
rz(-0.51760537) q[3];
sx q[3];
rz(2.8235265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.033325) q[2];
sx q[2];
rz(-2.6259618) q[2];
sx q[2];
rz(-2.2132614) q[2];
rz(0.27836529) q[3];
sx q[3];
rz(-1.3185578) q[3];
sx q[3];
rz(-0.068647169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6081029) q[0];
sx q[0];
rz(-2.6187496) q[0];
sx q[0];
rz(-2.9485517) q[0];
rz(-0.67735425) q[1];
sx q[1];
rz(-0.028379863) q[1];
sx q[1];
rz(1.631564) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91514912) q[0];
sx q[0];
rz(-2.9993176) q[0];
sx q[0];
rz(-0.14173843) q[0];
rz(-2.0329934) q[2];
sx q[2];
rz(-1.6110366) q[2];
sx q[2];
rz(1.6176335) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40367605) q[1];
sx q[1];
rz(-2.3737646) q[1];
sx q[1];
rz(1.9270792) q[1];
rz(2.5872748) q[3];
sx q[3];
rz(-0.80762562) q[3];
sx q[3];
rz(-0.19632158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7789) q[2];
sx q[2];
rz(-2.429246) q[2];
sx q[2];
rz(-0.98711291) q[2];
rz(-1.0846064) q[3];
sx q[3];
rz(-2.5955213) q[3];
sx q[3];
rz(-2.5206821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5622332) q[0];
sx q[0];
rz(-2.1491304) q[0];
sx q[0];
rz(1.5935422) q[0];
rz(0.69001895) q[1];
sx q[1];
rz(-3.1179805) q[1];
sx q[1];
rz(-1.4366368) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59047506) q[0];
sx q[0];
rz(-0.52177934) q[0];
sx q[0];
rz(0.087561122) q[0];
rz(-pi) q[1];
rz(1.5927926) q[2];
sx q[2];
rz(-2.2554033) q[2];
sx q[2];
rz(1.2134827) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1664983) q[1];
sx q[1];
rz(-1.0294895) q[1];
sx q[1];
rz(0.06947738) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2728746) q[3];
sx q[3];
rz(-2.0092589) q[3];
sx q[3];
rz(2.7496453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7684266) q[2];
sx q[2];
rz(-1.2370141) q[2];
sx q[2];
rz(0.58488673) q[2];
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
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5698513) q[0];
sx q[0];
rz(-2.2241156) q[0];
sx q[0];
rz(1.5745987) q[0];
rz(2.6245116) q[1];
sx q[1];
rz(-2.2061429) q[1];
sx q[1];
rz(1.2665952) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5735949) q[0];
sx q[0];
rz(-1.5657775) q[0];
sx q[0];
rz(-1.5742334) q[0];
rz(-2.1257519) q[2];
sx q[2];
rz(-2.2948143) q[2];
sx q[2];
rz(-0.99363771) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4609485) q[1];
sx q[1];
rz(-1.0310182) q[1];
sx q[1];
rz(-2.7311027) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4483767) q[3];
sx q[3];
rz(-1.4493692) q[3];
sx q[3];
rz(-1.9521015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7685585) q[2];
sx q[2];
rz(-2.3410102) q[2];
sx q[2];
rz(1.2726834) q[2];
rz(2.7019165) q[3];
sx q[3];
rz(-1.1594073) q[3];
sx q[3];
rz(-3.0769949) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8398297) q[0];
sx q[0];
rz(-0.01627144) q[0];
sx q[0];
rz(2.8453258) q[0];
rz(-0.06427327) q[1];
sx q[1];
rz(-1.4646894) q[1];
sx q[1];
rz(-1.6305264) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6600457) q[0];
sx q[0];
rz(-1.7007977) q[0];
sx q[0];
rz(1.2370716) q[0];
rz(0.6728386) q[2];
sx q[2];
rz(-1.7305859) q[2];
sx q[2];
rz(-0.21110134) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.293773) q[1];
sx q[1];
rz(-0.34983954) q[1];
sx q[1];
rz(-1.9560676) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6919208) q[3];
sx q[3];
rz(-1.559404) q[3];
sx q[3];
rz(-0.73339547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5955547) q[2];
sx q[2];
rz(-1.4348571) q[2];
sx q[2];
rz(1.1240553) q[2];
rz(1.837364) q[3];
sx q[3];
rz(-2.8458457) q[3];
sx q[3];
rz(-3.0834294) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27125636) q[0];
sx q[0];
rz(-2.3503292) q[0];
sx q[0];
rz(-1.9747718) q[0];
rz(1.5406436) q[1];
sx q[1];
rz(-0.29778844) q[1];
sx q[1];
rz(1.3265532) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1276774) q[0];
sx q[0];
rz(-1.8155037) q[0];
sx q[0];
rz(-0.6689594) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9341078) q[2];
sx q[2];
rz(-1.187154) q[2];
sx q[2];
rz(0.07019474) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.1804114) q[1];
sx q[1];
rz(-1.7827534) q[1];
sx q[1];
rz(-0.69976626) q[1];
rz(-1.972057) q[3];
sx q[3];
rz(-1.6607666) q[3];
sx q[3];
rz(-2.5800362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2549609) q[2];
sx q[2];
rz(-0.16356629) q[2];
sx q[2];
rz(-1.6939885) q[2];
rz(-0.12668954) q[3];
sx q[3];
rz(-1.4796175) q[3];
sx q[3];
rz(-2.0255585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3348677) q[0];
sx q[0];
rz(-1.3074449) q[0];
sx q[0];
rz(1.773651) q[0];
rz(1.5927636) q[1];
sx q[1];
rz(-2.3130885) q[1];
sx q[1];
rz(0.15514506) q[1];
rz(-2.6947179) q[2];
sx q[2];
rz(-1.4021654) q[2];
sx q[2];
rz(-1.5709102) q[2];
rz(2.1550989) q[3];
sx q[3];
rz(-1.7135847) q[3];
sx q[3];
rz(1.8691487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
