OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61385566) q[0];
sx q[0];
rz(4.6392415) q[0];
sx q[0];
rz(8.5949329) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(-0.87632626) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1286368) q[0];
sx q[0];
rz(-2.6814333) q[0];
sx q[0];
rz(-0.53928661) q[0];
x q[1];
rz(0.19199065) q[2];
sx q[2];
rz(-0.99388323) q[2];
sx q[2];
rz(2.752395) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10325601) q[1];
sx q[1];
rz(-1.5557319) q[1];
sx q[1];
rz(1.2407833) q[1];
rz(-pi) q[2];
x q[2];
rz(1.417744) q[3];
sx q[3];
rz(-1.1682086) q[3];
sx q[3];
rz(0.11868417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.17065389) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(-2.412964) q[2];
rz(0.5209926) q[3];
sx q[3];
rz(-2.1803768) q[3];
sx q[3];
rz(-2.9339824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8347297) q[0];
sx q[0];
rz(-1.9711718) q[0];
sx q[0];
rz(-1.1215425) q[0];
rz(-0.25575486) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(2.2671525) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5912936) q[0];
sx q[0];
rz(-1.8739788) q[0];
sx q[0];
rz(-1.1217872) q[0];
x q[1];
rz(-0.15408709) q[2];
sx q[2];
rz(-2.140897) q[2];
sx q[2];
rz(-2.4947583) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.16972152) q[1];
sx q[1];
rz(-1.1217146) q[1];
sx q[1];
rz(2.4596618) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3014684) q[3];
sx q[3];
rz(-1.519031) q[3];
sx q[3];
rz(-2.2477637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4008537) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(2.7056616) q[2];
rz(2.46051) q[3];
sx q[3];
rz(-0.77107945) q[3];
sx q[3];
rz(-2.7387103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23713672) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(-1.0748192) q[0];
rz(2.3020321) q[1];
sx q[1];
rz(-2.3222175) q[1];
sx q[1];
rz(0.39594617) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1186737) q[0];
sx q[0];
rz(-0.6624822) q[0];
sx q[0];
rz(0.89612095) q[0];
x q[1];
rz(-1.2767407) q[2];
sx q[2];
rz(-1.5573504) q[2];
sx q[2];
rz(-2.3384192) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6763941) q[1];
sx q[1];
rz(-0.68323831) q[1];
sx q[1];
rz(3.1264683) q[1];
x q[2];
rz(1.5474165) q[3];
sx q[3];
rz(-1.8262719) q[3];
sx q[3];
rz(-1.3083096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5376771) q[2];
sx q[2];
rz(-1.8751514) q[2];
sx q[2];
rz(-0.23920693) q[2];
rz(-3.0662597) q[3];
sx q[3];
rz(-2.0276666) q[3];
sx q[3];
rz(-1.8384365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72162119) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(-2.3216632) q[0];
rz(-2.6539102) q[1];
sx q[1];
rz(-2.2380424) q[1];
sx q[1];
rz(0.23342361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.061302) q[0];
sx q[0];
rz(-1.4740605) q[0];
sx q[0];
rz(-3.0808099) q[0];
rz(-pi) q[1];
rz(0.31077023) q[2];
sx q[2];
rz(-0.66590532) q[2];
sx q[2];
rz(-0.13194612) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1366795) q[1];
sx q[1];
rz(-2.3466952) q[1];
sx q[1];
rz(-2.6585048) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38846429) q[3];
sx q[3];
rz(-2.3187175) q[3];
sx q[3];
rz(-0.32271656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0908115) q[2];
sx q[2];
rz(-1.8596785) q[2];
sx q[2];
rz(-2.3941669) q[2];
rz(-0.22339544) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(-0.24211611) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6915879) q[0];
sx q[0];
rz(-1.7084028) q[0];
sx q[0];
rz(0.064237021) q[0];
rz(-0.94379395) q[1];
sx q[1];
rz(-2.4143024) q[1];
sx q[1];
rz(2.3805526) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9984263) q[0];
sx q[0];
rz(-2.879062) q[0];
sx q[0];
rz(-1.4545928) q[0];
rz(1.7663076) q[2];
sx q[2];
rz(-1.7695145) q[2];
sx q[2];
rz(0.18713258) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4796175) q[1];
sx q[1];
rz(-1.8912589) q[1];
sx q[1];
rz(-1.596405) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8588572) q[3];
sx q[3];
rz(-1.8394107) q[3];
sx q[3];
rz(2.3408567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3349907) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(-0.09207329) q[2];
rz(-0.66172415) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(1.3823284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(0.46835607) q[0];
sx q[0];
rz(-1.8259003) q[0];
sx q[0];
rz(-3.1150505) q[0];
rz(2.2684855) q[1];
sx q[1];
rz(-2.006242) q[1];
sx q[1];
rz(0.32593265) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78544261) q[0];
sx q[0];
rz(-2.4768562) q[0];
sx q[0];
rz(-0.63389969) q[0];
rz(-pi) q[1];
rz(-1.2665777) q[2];
sx q[2];
rz(-1.1433257) q[2];
sx q[2];
rz(1.3104591) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0210452) q[1];
sx q[1];
rz(-1.7056587) q[1];
sx q[1];
rz(-0.30121505) q[1];
x q[2];
rz(-3.0405007) q[3];
sx q[3];
rz(-1.5540082) q[3];
sx q[3];
rz(1.8198131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.59763336) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(0.65417543) q[2];
rz(-1.4298965) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(-0.54106075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27286801) q[0];
sx q[0];
rz(-1.6690212) q[0];
sx q[0];
rz(-2.4196999) q[0];
rz(1.4121274) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(-0.11925764) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3463466) q[0];
sx q[0];
rz(-0.7659142) q[0];
sx q[0];
rz(-2.7555097) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7330405) q[2];
sx q[2];
rz(-1.679323) q[2];
sx q[2];
rz(-0.42019368) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6357928) q[1];
sx q[1];
rz(-2.0675106) q[1];
sx q[1];
rz(-2.990681) q[1];
x q[2];
rz(-2.7151832) q[3];
sx q[3];
rz(-0.75424131) q[3];
sx q[3];
rz(-2.756556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44935903) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(1.3593486) q[2];
rz(2.3826777) q[3];
sx q[3];
rz(-2.9000498) q[3];
sx q[3];
rz(-0.53708491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3100202) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(2.0781562) q[0];
rz(-2.8670782) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(-2.2559821) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5078686) q[0];
sx q[0];
rz(-1.3390203) q[0];
sx q[0];
rz(-0.50110441) q[0];
rz(1.5947072) q[2];
sx q[2];
rz(-0.57966053) q[2];
sx q[2];
rz(-1.5469345) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9988261) q[1];
sx q[1];
rz(-2.4803659) q[1];
sx q[1];
rz(-0.89894526) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3732244) q[3];
sx q[3];
rz(-2.1185015) q[3];
sx q[3];
rz(2.2636569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4132335) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(-2.0098861) q[2];
rz(1.0845832) q[3];
sx q[3];
rz(-1.0794493) q[3];
sx q[3];
rz(-1.2148946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.8383012) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(1.0820748) q[0];
rz(1.2754296) q[1];
sx q[1];
rz(-1.0042896) q[1];
sx q[1];
rz(-1.1358322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1866859) q[0];
sx q[0];
rz(-1.0404772) q[0];
sx q[0];
rz(0.19219877) q[0];
rz(-2.7068107) q[2];
sx q[2];
rz(-1.5506622) q[2];
sx q[2];
rz(2.9895363) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75479924) q[1];
sx q[1];
rz(-2.7126185) q[1];
sx q[1];
rz(-3.0241443) q[1];
rz(-pi) q[2];
rz(-0.42580749) q[3];
sx q[3];
rz(-0.96447456) q[3];
sx q[3];
rz(-3.1283034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11848005) q[2];
sx q[2];
rz(-1.8818405) q[2];
sx q[2];
rz(-2.2686968) q[2];
rz(-0.84351271) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(-0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2492367) q[0];
sx q[0];
rz(-1.7974239) q[0];
sx q[0];
rz(-1.2783485) q[0];
rz(2.1168013) q[1];
sx q[1];
rz(-2.0147851) q[1];
sx q[1];
rz(1.1970253) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7322757) q[0];
sx q[0];
rz(-1.1336375) q[0];
sx q[0];
rz(0.93426312) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.808666) q[2];
sx q[2];
rz(-1.4537721) q[2];
sx q[2];
rz(-0.34840096) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5572284) q[1];
sx q[1];
rz(-0.96712501) q[1];
sx q[1];
rz(-0.61324688) q[1];
x q[2];
rz(-0.50935575) q[3];
sx q[3];
rz(-1.2378113) q[3];
sx q[3];
rz(-1.8978564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3060351) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(0.79997921) q[2];
rz(-1.9647313) q[3];
sx q[3];
rz(-1.4326982) q[3];
sx q[3];
rz(-2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70893127) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(2.6196383) q[1];
sx q[1];
rz(-2.3028761) q[1];
sx q[1];
rz(0.16470673) q[1];
rz(1.21571) q[2];
sx q[2];
rz(-1.9234895) q[2];
sx q[2];
rz(-1.1870155) q[2];
rz(2.5064777) q[3];
sx q[3];
rz(-1.0233581) q[3];
sx q[3];
rz(-0.81627853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
