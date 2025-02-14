OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47867632) q[0];
sx q[0];
rz(-1.588568) q[0];
sx q[0];
rz(1.8832062) q[0];
rz(-5.0212669) q[1];
sx q[1];
rz(6.8016383) q[1];
sx q[1];
rz(6.3432884) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4661048) q[0];
sx q[0];
rz(-1.4805313) q[0];
sx q[0];
rz(1.5399) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8608039) q[2];
sx q[2];
rz(-0.5092237) q[2];
sx q[2];
rz(-1.4851242) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2025238) q[1];
sx q[1];
rz(-1.2092522) q[1];
sx q[1];
rz(0.099379813) q[1];
rz(-pi) q[2];
rz(-2.2047049) q[3];
sx q[3];
rz(-0.70073116) q[3];
sx q[3];
rz(-2.7752293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1897159) q[2];
sx q[2];
rz(-0.78236255) q[2];
sx q[2];
rz(-0.18048364) q[2];
rz(0.30432025) q[3];
sx q[3];
rz(-2.2173827) q[3];
sx q[3];
rz(2.9144104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75905281) q[0];
sx q[0];
rz(-0.57342425) q[0];
sx q[0];
rz(-0.10391129) q[0];
rz(2.3567764) q[1];
sx q[1];
rz(-1.3094614) q[1];
sx q[1];
rz(-0.58194247) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3481929) q[0];
sx q[0];
rz(-1.5485067) q[0];
sx q[0];
rz(2.6096056) q[0];
rz(-pi) q[1];
rz(-2.632276) q[2];
sx q[2];
rz(-1.881372) q[2];
sx q[2];
rz(0.79307014) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4247893) q[1];
sx q[1];
rz(-1.5967073) q[1];
sx q[1];
rz(1.0938563) q[1];
rz(0.76762565) q[3];
sx q[3];
rz(-0.54477845) q[3];
sx q[3];
rz(1.2380074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4265784) q[2];
sx q[2];
rz(-2.7535186) q[2];
sx q[2];
rz(2.1132054) q[2];
rz(-0.55108023) q[3];
sx q[3];
rz(-1.9469399) q[3];
sx q[3];
rz(0.50857956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5904163) q[0];
sx q[0];
rz(-1.6352147) q[0];
sx q[0];
rz(-1.2580385) q[0];
rz(-2.538077) q[1];
sx q[1];
rz(-1.3110833) q[1];
sx q[1];
rz(-2.9579732) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4371673) q[0];
sx q[0];
rz(-1.1609623) q[0];
sx q[0];
rz(-2.8528575) q[0];
rz(2.6607291) q[2];
sx q[2];
rz(-1.4521952) q[2];
sx q[2];
rz(3.1080217) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.10079846) q[1];
sx q[1];
rz(-1.14193) q[1];
sx q[1];
rz(2.3593192) q[1];
rz(-1.612099) q[3];
sx q[3];
rz(-1.777711) q[3];
sx q[3];
rz(-2.4504091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.52643481) q[2];
sx q[2];
rz(-1.7650471) q[2];
sx q[2];
rz(-0.60205013) q[2];
rz(2.708882) q[3];
sx q[3];
rz(-1.5940462) q[3];
sx q[3];
rz(-1.4759147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7567619) q[0];
sx q[0];
rz(-0.47465208) q[0];
sx q[0];
rz(-0.62622825) q[0];
rz(-1.2681883) q[1];
sx q[1];
rz(-1.405193) q[1];
sx q[1];
rz(-0.38571206) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.884366) q[0];
sx q[0];
rz(-2.2472839) q[0];
sx q[0];
rz(-2.1348597) q[0];
rz(-pi) q[1];
rz(1.1168295) q[2];
sx q[2];
rz(-2.1056643) q[2];
sx q[2];
rz(1.7810389) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0929132) q[1];
sx q[1];
rz(-1.6374705) q[1];
sx q[1];
rz(-2.9408022) q[1];
rz(-pi) q[2];
rz(-1.0276876) q[3];
sx q[3];
rz(-1.6754136) q[3];
sx q[3];
rz(-1.0868768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4422153) q[2];
sx q[2];
rz(-2.0065887) q[2];
sx q[2];
rz(0.29083148) q[2];
rz(1.1902827) q[3];
sx q[3];
rz(-0.61605993) q[3];
sx q[3];
rz(-2.0975838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46665835) q[0];
sx q[0];
rz(-0.062677296) q[0];
sx q[0];
rz(-1.689893) q[0];
rz(-2.2845204) q[1];
sx q[1];
rz(-1.2985726) q[1];
sx q[1];
rz(2.7136386) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65954094) q[0];
sx q[0];
rz(-2.6845346) q[0];
sx q[0];
rz(1.8607451) q[0];
rz(-pi) q[1];
rz(2.4026947) q[2];
sx q[2];
rz(-1.1797) q[2];
sx q[2];
rz(-1.0223688) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2460096) q[1];
sx q[1];
rz(-1.2162788) q[1];
sx q[1];
rz(-2.406723) q[1];
rz(0.69231838) q[3];
sx q[3];
rz(-2.025617) q[3];
sx q[3];
rz(-2.2934283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2032623) q[2];
sx q[2];
rz(-2.0945956) q[2];
sx q[2];
rz(-0.14673512) q[2];
rz(-2.9444368) q[3];
sx q[3];
rz(-1.6517755) q[3];
sx q[3];
rz(2.3728235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7305304) q[0];
sx q[0];
rz(-1.0163607) q[0];
sx q[0];
rz(2.9732669) q[0];
rz(0.39237818) q[1];
sx q[1];
rz(-0.43586755) q[1];
sx q[1];
rz(1.6900774) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4844947) q[0];
sx q[0];
rz(-1.5494073) q[0];
sx q[0];
rz(1.1339897) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85003995) q[2];
sx q[2];
rz(-1.7687359) q[2];
sx q[2];
rz(-2.4356761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2830303) q[1];
sx q[1];
rz(-1.4301738) q[1];
sx q[1];
rz(2.3909516) q[1];
rz(1.1216954) q[3];
sx q[3];
rz(-0.83415745) q[3];
sx q[3];
rz(0.68461217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.77439857) q[2];
sx q[2];
rz(-2.6937679) q[2];
sx q[2];
rz(1.1773342) q[2];
rz(0.92159739) q[3];
sx q[3];
rz(-0.84601837) q[3];
sx q[3];
rz(0.53068501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.069799066) q[0];
sx q[0];
rz(-1.256218) q[0];
sx q[0];
rz(-1.0569093) q[0];
rz(2.9131556) q[1];
sx q[1];
rz(-2.3616796) q[1];
sx q[1];
rz(2.8905919) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91755501) q[0];
sx q[0];
rz(-0.38781181) q[0];
sx q[0];
rz(-0.73407956) q[0];
rz(-pi) q[1];
rz(0.86871882) q[2];
sx q[2];
rz(-1.0660604) q[2];
sx q[2];
rz(0.87930337) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2072548) q[1];
sx q[1];
rz(-1.435507) q[1];
sx q[1];
rz(0.16442085) q[1];
rz(-0.51336175) q[3];
sx q[3];
rz(-1.2436449) q[3];
sx q[3];
rz(1.038643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2061657) q[2];
sx q[2];
rz(-0.82038227) q[2];
sx q[2];
rz(2.7247735) q[2];
rz(-0.56944141) q[3];
sx q[3];
rz(-2.0408401) q[3];
sx q[3];
rz(-2.4429564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57132974) q[0];
sx q[0];
rz(-1.2705734) q[0];
sx q[0];
rz(-2.6369693) q[0];
rz(0.2568256) q[1];
sx q[1];
rz(-0.80373126) q[1];
sx q[1];
rz(1.1508734) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68783142) q[0];
sx q[0];
rz(-2.6021299) q[0];
sx q[0];
rz(-1.1629172) q[0];
rz(-pi) q[1];
rz(2.2058401) q[2];
sx q[2];
rz(-1.2798736) q[2];
sx q[2];
rz(2.3036602) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2765892) q[1];
sx q[1];
rz(-2.114608) q[1];
sx q[1];
rz(1.0629884) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7959149) q[3];
sx q[3];
rz(-0.44000235) q[3];
sx q[3];
rz(-1.0087165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5095832) q[2];
sx q[2];
rz(-2.7183967) q[2];
sx q[2];
rz(-0.43935856) q[2];
rz(-1.1257233) q[3];
sx q[3];
rz(-1.537354) q[3];
sx q[3];
rz(-1.2996947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87089649) q[0];
sx q[0];
rz(-0.82042158) q[0];
sx q[0];
rz(2.6743555) q[0];
rz(2.1029419) q[1];
sx q[1];
rz(-0.50913441) q[1];
sx q[1];
rz(1.6516623) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91083807) q[0];
sx q[0];
rz(-1.0808792) q[0];
sx q[0];
rz(1.7696437) q[0];
x q[1];
rz(2.5697904) q[2];
sx q[2];
rz(-0.17388136) q[2];
sx q[2];
rz(1.7553365) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4626325) q[1];
sx q[1];
rz(-2.3255146) q[1];
sx q[1];
rz(-0.036463375) q[1];
x q[2];
rz(-2.49754) q[3];
sx q[3];
rz(-0.67094147) q[3];
sx q[3];
rz(1.3356127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.44522875) q[2];
sx q[2];
rz(-2.7859521) q[2];
sx q[2];
rz(1.5544372) q[2];
rz(0.98443952) q[3];
sx q[3];
rz(-0.90960228) q[3];
sx q[3];
rz(-1.4136774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54870355) q[0];
sx q[0];
rz(-0.84243542) q[0];
sx q[0];
rz(-0.62826759) q[0];
rz(-2.9382622) q[1];
sx q[1];
rz(-1.1140099) q[1];
sx q[1];
rz(2.5471953) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0258163) q[0];
sx q[0];
rz(-2.5133142) q[0];
sx q[0];
rz(1.9253057) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6028676) q[2];
sx q[2];
rz(-1.3817412) q[2];
sx q[2];
rz(-2.3843887) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.17972936) q[1];
sx q[1];
rz(-0.41060251) q[1];
sx q[1];
rz(-0.73335464) q[1];
x q[2];
rz(1.6420308) q[3];
sx q[3];
rz(-2.1826577) q[3];
sx q[3];
rz(2.9737986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2079042) q[2];
sx q[2];
rz(-2.4173357) q[2];
sx q[2];
rz(0.17624632) q[2];
rz(1.8267953) q[3];
sx q[3];
rz(-0.81614554) q[3];
sx q[3];
rz(2.3740681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15277302) q[0];
sx q[0];
rz(-1.7012699) q[0];
sx q[0];
rz(-1.758601) q[0];
rz(3.026961) q[1];
sx q[1];
rz(-0.75420598) q[1];
sx q[1];
rz(2.4892714) q[1];
rz(0.23723142) q[2];
sx q[2];
rz(-1.5330412) q[2];
sx q[2];
rz(1.588149) q[2];
rz(-0.38269855) q[3];
sx q[3];
rz(-2.0321587) q[3];
sx q[3];
rz(1.7176499) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
