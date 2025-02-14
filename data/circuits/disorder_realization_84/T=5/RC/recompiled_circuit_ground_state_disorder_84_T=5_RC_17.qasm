OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.428838) q[0];
sx q[0];
rz(-1.5313671) q[0];
sx q[0];
rz(2.03696) q[0];
rz(1.334335) q[1];
sx q[1];
rz(-2.4185138) q[1];
sx q[1];
rz(1.877797) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1570779) q[0];
sx q[0];
rz(-2.4369168) q[0];
sx q[0];
rz(-1.1873755) q[0];
x q[1];
rz(-0.96226389) q[2];
sx q[2];
rz(-2.4748286) q[2];
sx q[2];
rz(-2.5279669) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7967148) q[1];
sx q[1];
rz(-1.2462707) q[1];
sx q[1];
rz(-0.085436324) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8975211) q[3];
sx q[3];
rz(-1.1435978) q[3];
sx q[3];
rz(1.106316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9870712) q[2];
sx q[2];
rz(-2.0646586) q[2];
sx q[2];
rz(1.3686352) q[2];
rz(0.23855071) q[3];
sx q[3];
rz(-0.41540256) q[3];
sx q[3];
rz(0.11428782) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7411165) q[0];
sx q[0];
rz(-1.1392765) q[0];
sx q[0];
rz(1.9912632) q[0];
rz(-0.087609619) q[1];
sx q[1];
rz(-1.8116415) q[1];
sx q[1];
rz(1.5709343) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51975361) q[0];
sx q[0];
rz(-1.4645508) q[0];
sx q[0];
rz(2.6563717) q[0];
rz(1.6341799) q[2];
sx q[2];
rz(-1.7492314) q[2];
sx q[2];
rz(-3.1270535) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0391866) q[1];
sx q[1];
rz(-0.79823179) q[1];
sx q[1];
rz(-2.185553) q[1];
x q[2];
rz(2.0474954) q[3];
sx q[3];
rz(-2.3040651) q[3];
sx q[3];
rz(0.9707211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34438434) q[2];
sx q[2];
rz(-0.71343652) q[2];
sx q[2];
rz(0.1864645) q[2];
rz(-2.3421085) q[3];
sx q[3];
rz(-1.5954285) q[3];
sx q[3];
rz(-0.98108393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28476533) q[0];
sx q[0];
rz(-1.7600049) q[0];
sx q[0];
rz(-3.0322266) q[0];
rz(-2.5378387) q[1];
sx q[1];
rz(-1.8021288) q[1];
sx q[1];
rz(1.0341136) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17149481) q[0];
sx q[0];
rz(-1.437828) q[0];
sx q[0];
rz(1.7823969) q[0];
rz(-1.500152) q[2];
sx q[2];
rz(-1.4817258) q[2];
sx q[2];
rz(-1.8235109) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.32419606) q[1];
sx q[1];
rz(-1.6764371) q[1];
sx q[1];
rz(-0.28075851) q[1];
rz(-pi) q[2];
rz(-1.4294741) q[3];
sx q[3];
rz(-0.80414786) q[3];
sx q[3];
rz(0.49403163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5321396) q[2];
sx q[2];
rz(-0.9708465) q[2];
sx q[2];
rz(-2.5965221) q[2];
rz(0.80101454) q[3];
sx q[3];
rz(-2.2478734) q[3];
sx q[3];
rz(3.0494704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48701778) q[0];
sx q[0];
rz(-1.6681404) q[0];
sx q[0];
rz(0.45904485) q[0];
rz(2.0252939) q[1];
sx q[1];
rz(-2.3479159) q[1];
sx q[1];
rz(2.3150516) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61454489) q[0];
sx q[0];
rz(-1.2189208) q[0];
sx q[0];
rz(-1.1220758) q[0];
rz(0.79016308) q[2];
sx q[2];
rz(-0.1383257) q[2];
sx q[2];
rz(-2.9776855) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.53163995) q[1];
sx q[1];
rz(-2.5594257) q[1];
sx q[1];
rz(-1.9332631) q[1];
rz(-pi) q[2];
rz(1.4524121) q[3];
sx q[3];
rz(-1.9031798) q[3];
sx q[3];
rz(0.82086241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7894342) q[2];
sx q[2];
rz(-0.23098478) q[2];
sx q[2];
rz(-0.74971548) q[2];
rz(1.1226783) q[3];
sx q[3];
rz(-1.2238945) q[3];
sx q[3];
rz(-0.96021715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.66626755) q[0];
sx q[0];
rz(-0.98639494) q[0];
sx q[0];
rz(3.0295897) q[0];
rz(2.9169967) q[1];
sx q[1];
rz(-1.1677531) q[1];
sx q[1];
rz(-2.3602643) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8664426) q[0];
sx q[0];
rz(-2.7366182) q[0];
sx q[0];
rz(-1.6385796) q[0];
rz(-2.8728169) q[2];
sx q[2];
rz(-2.7999863) q[2];
sx q[2];
rz(-2.9363321) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.71415662) q[1];
sx q[1];
rz(-0.86227741) q[1];
sx q[1];
rz(-0.57517131) q[1];
x q[2];
rz(2.1634401) q[3];
sx q[3];
rz(-1.4362122) q[3];
sx q[3];
rz(2.5434567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8319228) q[2];
sx q[2];
rz(-0.90193844) q[2];
sx q[2];
rz(0.93264467) q[2];
rz(2.0617088) q[3];
sx q[3];
rz(-2.6974758) q[3];
sx q[3];
rz(-3.1058969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90750736) q[0];
sx q[0];
rz(-2.0103173) q[0];
sx q[0];
rz(-1.3908516) q[0];
rz(2.8098409) q[1];
sx q[1];
rz(-1.2650047) q[1];
sx q[1];
rz(1.5914241) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0557779) q[0];
sx q[0];
rz(-1.6059438) q[0];
sx q[0];
rz(-2.2136392) q[0];
x q[1];
rz(1.5339025) q[2];
sx q[2];
rz(-1.1114745) q[2];
sx q[2];
rz(1.5609891) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5065803) q[1];
sx q[1];
rz(-1.6939347) q[1];
sx q[1];
rz(-1.413373) q[1];
rz(-1.9523296) q[3];
sx q[3];
rz(-1.23151) q[3];
sx q[3];
rz(1.3374223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8288237) q[2];
sx q[2];
rz(-0.88461107) q[2];
sx q[2];
rz(-1.5213607) q[2];
rz(2.17365) q[3];
sx q[3];
rz(-1.5827551) q[3];
sx q[3];
rz(0.91606417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62436002) q[0];
sx q[0];
rz(-0.42651287) q[0];
sx q[0];
rz(0.17159167) q[0];
rz(2.102237) q[1];
sx q[1];
rz(-1.8828705) q[1];
sx q[1];
rz(-1.7003869) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0943308) q[0];
sx q[0];
rz(-2.7956595) q[0];
sx q[0];
rz(0.37389116) q[0];
rz(-pi) q[1];
rz(-0.19564678) q[2];
sx q[2];
rz(-0.52426978) q[2];
sx q[2];
rz(1.153868) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4413195) q[1];
sx q[1];
rz(-2.0236349) q[1];
sx q[1];
rz(0.39703607) q[1];
x q[2];
rz(1.7923545) q[3];
sx q[3];
rz(-2.5668813) q[3];
sx q[3];
rz(1.2971085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.27199304) q[2];
sx q[2];
rz(-1.4078434) q[2];
sx q[2];
rz(-0.54455152) q[2];
rz(0.49992391) q[3];
sx q[3];
rz(-2.1130424) q[3];
sx q[3];
rz(0.47206363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2917824) q[0];
sx q[0];
rz(-0.53684679) q[0];
sx q[0];
rz(-2.612402) q[0];
rz(-2.5657907) q[1];
sx q[1];
rz(-1.2628097) q[1];
sx q[1];
rz(-1.1189438) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81572318) q[0];
sx q[0];
rz(-1.5191606) q[0];
sx q[0];
rz(-1.5918094) q[0];
rz(-pi) q[1];
rz(0.37092692) q[2];
sx q[2];
rz(-1.7417522) q[2];
sx q[2];
rz(-1.2990739) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0487489) q[1];
sx q[1];
rz(-2.1511937) q[1];
sx q[1];
rz(-1.012085) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0980074) q[3];
sx q[3];
rz(-1.7526748) q[3];
sx q[3];
rz(-1.393569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.240856) q[2];
sx q[2];
rz(-0.46515981) q[2];
sx q[2];
rz(-2.4294803) q[2];
rz(1.6093048) q[3];
sx q[3];
rz(-0.94523793) q[3];
sx q[3];
rz(-1.3473264) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3466472) q[0];
sx q[0];
rz(-2.0136588) q[0];
sx q[0];
rz(-2.0704863) q[0];
rz(1.2062997) q[1];
sx q[1];
rz(-2.1251528) q[1];
sx q[1];
rz(1.4869022) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7506055) q[0];
sx q[0];
rz(-1.3840527) q[0];
sx q[0];
rz(-2.4658527) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73545154) q[2];
sx q[2];
rz(-1.1142154) q[2];
sx q[2];
rz(0.18010715) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2675954) q[1];
sx q[1];
rz(-1.8615627) q[1];
sx q[1];
rz(-1.7520755) q[1];
rz(2.1110299) q[3];
sx q[3];
rz(-0.5117473) q[3];
sx q[3];
rz(-2.668021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6335166) q[2];
sx q[2];
rz(-0.8997007) q[2];
sx q[2];
rz(-3.0637975) q[2];
rz(2.9142694) q[3];
sx q[3];
rz(-1.1319755) q[3];
sx q[3];
rz(-1.0556861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8530497) q[0];
sx q[0];
rz(-0.74497861) q[0];
sx q[0];
rz(-0.37260923) q[0];
rz(0.44652069) q[1];
sx q[1];
rz(-1.4563072) q[1];
sx q[1];
rz(2.2850697) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4093709) q[0];
sx q[0];
rz(-1.6147095) q[0];
sx q[0];
rz(-0.056890566) q[0];
rz(-2.2394453) q[2];
sx q[2];
rz(-0.16996102) q[2];
sx q[2];
rz(-0.68357498) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0792993) q[1];
sx q[1];
rz(-2.7253299) q[1];
sx q[1];
rz(1.2757311) q[1];
rz(-pi) q[2];
rz(2.3395446) q[3];
sx q[3];
rz(-1.060876) q[3];
sx q[3];
rz(-2.0039441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38241688) q[2];
sx q[2];
rz(-2.0903812) q[2];
sx q[2];
rz(-0.55553931) q[2];
rz(2.7555452) q[3];
sx q[3];
rz(-1.4945533) q[3];
sx q[3];
rz(0.66421318) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1356708) q[0];
sx q[0];
rz(-1.6914524) q[0];
sx q[0];
rz(1.2941262) q[0];
rz(0.80264965) q[1];
sx q[1];
rz(-1.4603271) q[1];
sx q[1];
rz(-2.7535798) q[1];
rz(-1.6907202) q[2];
sx q[2];
rz(-1.220206) q[2];
sx q[2];
rz(-2.7239067) q[2];
rz(0.070128154) q[3];
sx q[3];
rz(-2.3002516) q[3];
sx q[3];
rz(-2.3723835) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
