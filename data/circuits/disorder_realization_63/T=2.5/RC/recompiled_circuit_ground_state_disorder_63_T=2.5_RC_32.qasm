OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11290057) q[0];
sx q[0];
rz(-0.24554645) q[0];
sx q[0];
rz(-0.052659642) q[0];
rz(1.117299) q[1];
sx q[1];
rz(-2.8007562) q[1];
sx q[1];
rz(-2.053082) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77086335) q[0];
sx q[0];
rz(-1.5077871) q[0];
sx q[0];
rz(1.9669238) q[0];
rz(-0.11012385) q[2];
sx q[2];
rz(-1.4297669) q[2];
sx q[2];
rz(1.5361496) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.12456914) q[1];
sx q[1];
rz(-0.93734159) q[1];
sx q[1];
rz(1.3658466) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9904706) q[3];
sx q[3];
rz(-1.5670781) q[3];
sx q[3];
rz(-2.9633455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5071621) q[2];
sx q[2];
rz(-1.7796702) q[2];
sx q[2];
rz(0.6081028) q[2];
rz(2.0173343) q[3];
sx q[3];
rz(-1.2048771) q[3];
sx q[3];
rz(0.89157909) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41481498) q[0];
sx q[0];
rz(-1.5730653) q[0];
sx q[0];
rz(-2.538105) q[0];
rz(-2.6314349) q[1];
sx q[1];
rz(-2.3699103) q[1];
sx q[1];
rz(-1.0268432) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3393766) q[0];
sx q[0];
rz(-0.23288865) q[0];
sx q[0];
rz(-2.4075137) q[0];
x q[1];
rz(0.6925236) q[2];
sx q[2];
rz(-2.1998458) q[2];
sx q[2];
rz(0.72450996) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.75637792) q[1];
sx q[1];
rz(-0.88895386) q[1];
sx q[1];
rz(-2.4663976) q[1];
rz(-0.92225109) q[3];
sx q[3];
rz(-0.77510683) q[3];
sx q[3];
rz(0.67476455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5828731) q[2];
sx q[2];
rz(-1.9931953) q[2];
sx q[2];
rz(-0.64644512) q[2];
rz(1.2785814) q[3];
sx q[3];
rz(-0.94978142) q[3];
sx q[3];
rz(-1.6897374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.932514) q[0];
sx q[0];
rz(-0.13020733) q[0];
sx q[0];
rz(1.5721488) q[0];
rz(2.7945844) q[1];
sx q[1];
rz(-1.4748814) q[1];
sx q[1];
rz(-2.2775547) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0542674) q[0];
sx q[0];
rz(-2.1313166) q[0];
sx q[0];
rz(1.3232687) q[0];
x q[1];
rz(-2.6367497) q[2];
sx q[2];
rz(-2.6528203) q[2];
sx q[2];
rz(-0.29112838) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.62926312) q[1];
sx q[1];
rz(-2.1844668) q[1];
sx q[1];
rz(-0.83028173) q[1];
x q[2];
rz(2.0967617) q[3];
sx q[3];
rz(-1.6623896) q[3];
sx q[3];
rz(-2.4909702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.049909264) q[2];
sx q[2];
rz(-0.20541643) q[2];
sx q[2];
rz(-1.6845711) q[2];
rz(2.125804) q[3];
sx q[3];
rz(-2.6088645) q[3];
sx q[3];
rz(0.31639019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-2.9006573) q[0];
sx q[0];
rz(-2.7641986) q[0];
sx q[0];
rz(1.5628016) q[0];
rz(-2.0511625) q[1];
sx q[1];
rz(-0.54160392) q[1];
sx q[1];
rz(-1.1515559) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82842302) q[0];
sx q[0];
rz(-2.1439664) q[0];
sx q[0];
rz(1.8445119) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3811023) q[2];
sx q[2];
rz(-1.2142688) q[2];
sx q[2];
rz(-0.52887756) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87428625) q[1];
sx q[1];
rz(-1.6055672) q[1];
sx q[1];
rz(2.094993) q[1];
rz(1.6107481) q[3];
sx q[3];
rz(-2.5343347) q[3];
sx q[3];
rz(1.1049529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2942723) q[2];
sx q[2];
rz(-2.5966094) q[2];
sx q[2];
rz(-1.4997743) q[2];
rz(0.13449399) q[3];
sx q[3];
rz(-1.8650863) q[3];
sx q[3];
rz(-0.78181481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0791557) q[0];
sx q[0];
rz(-0.9117313) q[0];
sx q[0];
rz(2.675918) q[0];
rz(1.8648719) q[1];
sx q[1];
rz(-1.0036422) q[1];
sx q[1];
rz(-1.0328971) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0807652) q[0];
sx q[0];
rz(-1.5603335) q[0];
sx q[0];
rz(1.7288789) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6200142) q[2];
sx q[2];
rz(-0.17159941) q[2];
sx q[2];
rz(2.624334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3001662) q[1];
sx q[1];
rz(-2.4304217) q[1];
sx q[1];
rz(0.39022846) q[1];
rz(-pi) q[2];
rz(2.1139268) q[3];
sx q[3];
rz(-1.2333721) q[3];
sx q[3];
rz(1.8130482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5501962) q[2];
sx q[2];
rz(-1.952221) q[2];
sx q[2];
rz(-0.93185321) q[2];
rz(-0.094680928) q[3];
sx q[3];
rz(-0.90016142) q[3];
sx q[3];
rz(-1.3100821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4513627) q[0];
sx q[0];
rz(-2.9424423) q[0];
sx q[0];
rz(3.0354101) q[0];
rz(-0.55446082) q[1];
sx q[1];
rz(-0.78840557) q[1];
sx q[1];
rz(-1.5415446) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43801895) q[0];
sx q[0];
rz(-0.39123639) q[0];
sx q[0];
rz(-3.0783669) q[0];
rz(0.27435209) q[2];
sx q[2];
rz(-1.128528) q[2];
sx q[2];
rz(-1.2081462) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.627085) q[1];
sx q[1];
rz(-0.70620757) q[1];
sx q[1];
rz(0.18893623) q[1];
rz(-0.24644773) q[3];
sx q[3];
rz(-1.5165899) q[3];
sx q[3];
rz(-1.1185631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23325486) q[2];
sx q[2];
rz(-2.4007863) q[2];
sx q[2];
rz(1.4976658) q[2];
rz(2.6805367) q[3];
sx q[3];
rz(-0.65270972) q[3];
sx q[3];
rz(1.8259995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4696362) q[0];
sx q[0];
rz(-2.3905583) q[0];
sx q[0];
rz(-1.0787971) q[0];
rz(-0.59393334) q[1];
sx q[1];
rz(-2.1682231) q[1];
sx q[1];
rz(2.914391) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57892913) q[0];
sx q[0];
rz(-1.4592071) q[0];
sx q[0];
rz(-1.4212378) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.960177) q[2];
sx q[2];
rz(-2.7837278) q[2];
sx q[2];
rz(0.030454446) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.46887423) q[1];
sx q[1];
rz(-2.1203845) q[1];
sx q[1];
rz(-0.36512111) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9428857) q[3];
sx q[3];
rz(-1.477302) q[3];
sx q[3];
rz(-1.0013415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6195153) q[2];
sx q[2];
rz(-2.7232813) q[2];
sx q[2];
rz(-2.412001) q[2];
rz(1.6884165) q[3];
sx q[3];
rz(-1.6875234) q[3];
sx q[3];
rz(-2.5280473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59298092) q[0];
sx q[0];
rz(-0.91658968) q[0];
sx q[0];
rz(2.7899637) q[0];
rz(0.31373203) q[1];
sx q[1];
rz(-1.2064563) q[1];
sx q[1];
rz(-1.3765913) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2257654) q[0];
sx q[0];
rz(-1.6179913) q[0];
sx q[0];
rz(-2.7412358) q[0];
x q[1];
rz(-1.7369164) q[2];
sx q[2];
rz(-1.4567647) q[2];
sx q[2];
rz(1.5064552) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.541988) q[1];
sx q[1];
rz(-1.7270589) q[1];
sx q[1];
rz(-2.5708052) q[1];
x q[2];
rz(1.9405211) q[3];
sx q[3];
rz(-1.4784357) q[3];
sx q[3];
rz(-0.28932387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5967963) q[2];
sx q[2];
rz(-2.4185541) q[2];
sx q[2];
rz(1.6205622) q[2];
rz(-0.14791402) q[3];
sx q[3];
rz(-1.3444129) q[3];
sx q[3];
rz(1.3676876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77383298) q[0];
sx q[0];
rz(-1.8748883) q[0];
sx q[0];
rz(1.2441147) q[0];
rz(-1.4338088) q[1];
sx q[1];
rz(-1.5807187) q[1];
sx q[1];
rz(0.22020766) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4996195) q[0];
sx q[0];
rz(-0.73055482) q[0];
sx q[0];
rz(-2.1424908) q[0];
x q[1];
rz(-0.16248361) q[2];
sx q[2];
rz(-1.1097317) q[2];
sx q[2];
rz(1.7562255) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23277321) q[1];
sx q[1];
rz(-2.4202883) q[1];
sx q[1];
rz(-2.2386293) q[1];
rz(2.3784786) q[3];
sx q[3];
rz(-3.0060351) q[3];
sx q[3];
rz(-0.28597304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1844909) q[2];
sx q[2];
rz(-1.1471006) q[2];
sx q[2];
rz(-1.5000878) q[2];
rz(-0.084065048) q[3];
sx q[3];
rz(-1.8819239) q[3];
sx q[3];
rz(3.0715004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67749196) q[0];
sx q[0];
rz(-2.3832432) q[0];
sx q[0];
rz(-0.41859928) q[0];
rz(2.1981926) q[1];
sx q[1];
rz(-1.489233) q[1];
sx q[1];
rz(2.8387866) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3269791) q[0];
sx q[0];
rz(-1.4340861) q[0];
sx q[0];
rz(1.4025406) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19089107) q[2];
sx q[2];
rz(-1.396469) q[2];
sx q[2];
rz(3.0534923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8119366) q[1];
sx q[1];
rz(-1.4881663) q[1];
sx q[1];
rz(-1.5272115) q[1];
x q[2];
rz(3.0251011) q[3];
sx q[3];
rz(-1.2283652) q[3];
sx q[3];
rz(0.064700944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6454978) q[2];
sx q[2];
rz(-0.87254137) q[2];
sx q[2];
rz(-2.3910451) q[2];
rz(1.4613072) q[3];
sx q[3];
rz(-2.3716726) q[3];
sx q[3];
rz(-2.6245978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1207598) q[0];
sx q[0];
rz(-1.5789541) q[0];
sx q[0];
rz(-1.5776237) q[0];
rz(-1.9137406) q[1];
sx q[1];
rz(-2.6422983) q[1];
sx q[1];
rz(-1.9849389) q[1];
rz(-0.56671178) q[2];
sx q[2];
rz(-1.4162345) q[2];
sx q[2];
rz(-0.92879374) q[2];
rz(1.8932992) q[3];
sx q[3];
rz(-1.396198) q[3];
sx q[3];
rz(-0.085343501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
