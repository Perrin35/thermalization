OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6629163) q[0];
sx q[0];
rz(-1.5530246) q[0];
sx q[0];
rz(-1.8832062) q[0];
rz(1.2619184) q[1];
sx q[1];
rz(-2.6231397) q[1];
sx q[1];
rz(3.0814896) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10190554) q[0];
sx q[0];
rz(-1.6015669) q[0];
sx q[0];
rz(3.0512848) q[0];
x q[1];
rz(2.6491935) q[2];
sx q[2];
rz(-1.4352893) q[2];
sx q[2];
rz(-2.8092334) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5451193) q[1];
sx q[1];
rz(-1.6637322) q[1];
sx q[1];
rz(-1.93398) q[1];
x q[2];
rz(0.93688776) q[3];
sx q[3];
rz(-0.70073116) q[3];
sx q[3];
rz(-2.7752293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1897159) q[2];
sx q[2];
rz(-0.78236255) q[2];
sx q[2];
rz(0.18048364) q[2];
rz(0.30432025) q[3];
sx q[3];
rz(-0.92420998) q[3];
sx q[3];
rz(-2.9144104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75905281) q[0];
sx q[0];
rz(-2.5681684) q[0];
sx q[0];
rz(0.10391129) q[0];
rz(-0.78481627) q[1];
sx q[1];
rz(-1.8321313) q[1];
sx q[1];
rz(-2.5596502) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4020445) q[0];
sx q[0];
rz(-2.6091837) q[0];
sx q[0];
rz(-3.0976712) q[0];
rz(-pi) q[1];
rz(0.50931661) q[2];
sx q[2];
rz(-1.2602206) q[2];
sx q[2];
rz(-0.79307014) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.71680333) q[1];
sx q[1];
rz(-1.5967073) q[1];
sx q[1];
rz(-2.0477363) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76762565) q[3];
sx q[3];
rz(-0.54477845) q[3];
sx q[3];
rz(1.9035853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71501422) q[2];
sx q[2];
rz(-2.7535186) q[2];
sx q[2];
rz(1.0283872) q[2];
rz(0.55108023) q[3];
sx q[3];
rz(-1.1946528) q[3];
sx q[3];
rz(0.50857956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5904163) q[0];
sx q[0];
rz(-1.506378) q[0];
sx q[0];
rz(1.8835541) q[0];
rz(-2.538077) q[1];
sx q[1];
rz(-1.8305093) q[1];
sx q[1];
rz(2.9579732) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98417898) q[0];
sx q[0];
rz(-1.3065728) q[0];
sx q[0];
rz(-1.1452894) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8894561) q[2];
sx q[2];
rz(-0.49415961) q[2];
sx q[2];
rz(-1.8273938) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0634707) q[1];
sx q[1];
rz(-2.2666711) q[1];
sx q[1];
rz(-2.1433926) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20708624) q[3];
sx q[3];
rz(-1.5303751) q[3];
sx q[3];
rz(-2.25349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6151578) q[2];
sx q[2];
rz(-1.7650471) q[2];
sx q[2];
rz(2.5395425) q[2];
rz(-2.708882) q[3];
sx q[3];
rz(-1.5940462) q[3];
sx q[3];
rz(-1.6656779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3848307) q[0];
sx q[0];
rz(-0.47465208) q[0];
sx q[0];
rz(2.5153644) q[0];
rz(1.8734044) q[1];
sx q[1];
rz(-1.405193) q[1];
sx q[1];
rz(2.7558806) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69067467) q[0];
sx q[0];
rz(-2.0008149) q[0];
sx q[0];
rz(0.7597835) q[0];
rz(1.1168295) q[2];
sx q[2];
rz(-1.0359284) q[2];
sx q[2];
rz(-1.7810389) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9358237) q[1];
sx q[1];
rz(-0.21142928) q[1];
sx q[1];
rz(0.32306674) q[1];
x q[2];
rz(2.1139051) q[3];
sx q[3];
rz(-1.6754136) q[3];
sx q[3];
rz(2.0547158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6993774) q[2];
sx q[2];
rz(-1.135004) q[2];
sx q[2];
rz(-0.29083148) q[2];
rz(-1.95131) q[3];
sx q[3];
rz(-0.61605993) q[3];
sx q[3];
rz(-2.0975838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6749343) q[0];
sx q[0];
rz(-0.062677296) q[0];
sx q[0];
rz(-1.4516996) q[0];
rz(-0.85707227) q[1];
sx q[1];
rz(-1.2985726) q[1];
sx q[1];
rz(-2.7136386) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1728522) q[0];
sx q[0];
rz(-1.6973054) q[0];
sx q[0];
rz(2.0111903) q[0];
x q[1];
rz(0.73889795) q[2];
sx q[2];
rz(-1.1797) q[2];
sx q[2];
rz(-2.1192239) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.9792069) q[1];
sx q[1];
rz(-0.89079327) q[1];
sx q[1];
rz(-1.1080145) q[1];
rz(-0.6537207) q[3];
sx q[3];
rz(-2.3343784) q[3];
sx q[3];
rz(1.2098055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9383303) q[2];
sx q[2];
rz(-1.0469971) q[2];
sx q[2];
rz(0.14673512) q[2];
rz(-0.19715582) q[3];
sx q[3];
rz(-1.6517755) q[3];
sx q[3];
rz(-2.3728235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41106221) q[0];
sx q[0];
rz(-2.125232) q[0];
sx q[0];
rz(-2.9732669) q[0];
rz(-0.39237818) q[1];
sx q[1];
rz(-0.43586755) q[1];
sx q[1];
rz(-1.6900774) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0453059) q[0];
sx q[0];
rz(-2.0074962) q[0];
sx q[0];
rz(-0.023604579) q[0];
rz(-1.8658357) q[2];
sx q[2];
rz(-0.74271281) q[2];
sx q[2];
rz(-2.4969522) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85856232) q[1];
sx q[1];
rz(-1.4301738) q[1];
sx q[1];
rz(-2.3909516) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78877641) q[3];
sx q[3];
rz(-1.8982072) q[3];
sx q[3];
rz(0.57306266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3671941) q[2];
sx q[2];
rz(-2.6937679) q[2];
sx q[2];
rz(1.9642584) q[2];
rz(2.2199953) q[3];
sx q[3];
rz(-2.2955743) q[3];
sx q[3];
rz(0.53068501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069799066) q[0];
sx q[0];
rz(-1.256218) q[0];
sx q[0];
rz(-2.0846833) q[0];
rz(-0.22843703) q[1];
sx q[1];
rz(-2.3616796) q[1];
sx q[1];
rz(2.8905919) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0989444) q[0];
sx q[0];
rz(-1.8269208) q[0];
sx q[0];
rz(2.8471208) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2784581) q[2];
sx q[2];
rz(-2.302711) q[2];
sx q[2];
rz(1.9307435) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9343379) q[1];
sx q[1];
rz(-1.435507) q[1];
sx q[1];
rz(-0.16442085) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5369473) q[3];
sx q[3];
rz(-2.5408158) q[3];
sx q[3];
rz(1.050211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.935427) q[2];
sx q[2];
rz(-0.82038227) q[2];
sx q[2];
rz(2.7247735) q[2];
rz(-2.5721512) q[3];
sx q[3];
rz(-1.1007525) q[3];
sx q[3];
rz(-2.4429564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57132974) q[0];
sx q[0];
rz(-1.8710192) q[0];
sx q[0];
rz(-2.6369693) q[0];
rz(-2.8847671) q[1];
sx q[1];
rz(-0.80373126) q[1];
sx q[1];
rz(1.1508734) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68783142) q[0];
sx q[0];
rz(-0.53946278) q[0];
sx q[0];
rz(1.1629172) q[0];
rz(-pi) q[1];
rz(0.93575259) q[2];
sx q[2];
rz(-1.861719) q[2];
sx q[2];
rz(-0.83793241) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.013864036) q[1];
sx q[1];
rz(-1.999966) q[1];
sx q[1];
rz(0.60529373) q[1];
rz(-pi) q[2];
rz(-2.7959149) q[3];
sx q[3];
rz(-2.7015903) q[3];
sx q[3];
rz(-2.1328762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5095832) q[2];
sx q[2];
rz(-0.42319599) q[2];
sx q[2];
rz(0.43935856) q[2];
rz(-2.0158694) q[3];
sx q[3];
rz(-1.537354) q[3];
sx q[3];
rz(1.2996947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2706962) q[0];
sx q[0];
rz(-0.82042158) q[0];
sx q[0];
rz(0.46723715) q[0];
rz(-1.0386508) q[1];
sx q[1];
rz(-2.6324582) q[1];
sx q[1];
rz(-1.6516623) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3154513) q[0];
sx q[0];
rz(-2.6159161) q[0];
sx q[0];
rz(-2.7868411) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5697904) q[2];
sx q[2];
rz(-2.9677113) q[2];
sx q[2];
rz(1.3862561) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6789602) q[1];
sx q[1];
rz(-2.3255146) q[1];
sx q[1];
rz(-0.036463375) q[1];
x q[2];
rz(-2.0155679) q[3];
sx q[3];
rz(-1.0504587) q[3];
sx q[3];
rz(-2.570278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6963639) q[2];
sx q[2];
rz(-2.7859521) q[2];
sx q[2];
rz(-1.5544372) q[2];
rz(0.98443952) q[3];
sx q[3];
rz(-0.90960228) q[3];
sx q[3];
rz(1.7279153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54870355) q[0];
sx q[0];
rz(-2.2991572) q[0];
sx q[0];
rz(-0.62826759) q[0];
rz(0.20333044) q[1];
sx q[1];
rz(-2.0275828) q[1];
sx q[1];
rz(-2.5471953) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1157764) q[0];
sx q[0];
rz(-0.62827841) q[0];
sx q[0];
rz(1.9253057) q[0];
rz(-pi) q[1];
rz(-1.3514694) q[2];
sx q[2];
rz(-2.0989053) q[2];
sx q[2];
rz(-0.92544941) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.59692057) q[1];
sx q[1];
rz(-1.2697176) q[1];
sx q[1];
rz(1.8543509) q[1];
rz(0.61305586) q[3];
sx q[3];
rz(-1.5125015) q[3];
sx q[3];
rz(1.4439652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9336885) q[2];
sx q[2];
rz(-0.7242569) q[2];
sx q[2];
rz(0.17624632) q[2];
rz(1.8267953) q[3];
sx q[3];
rz(-2.3254471) q[3];
sx q[3];
rz(0.76752457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15277302) q[0];
sx q[0];
rz(-1.4403227) q[0];
sx q[0];
rz(1.3829917) q[0];
rz(-3.026961) q[1];
sx q[1];
rz(-2.3873867) q[1];
sx q[1];
rz(-0.65232123) q[1];
rz(-0.15936511) q[2];
sx q[2];
rz(-2.9014316) q[2];
sx q[2];
rz(0.17221774) q[2];
rz(0.38269855) q[3];
sx q[3];
rz(-1.1094339) q[3];
sx q[3];
rz(-1.4239428) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
