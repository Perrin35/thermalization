OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(2.7352754) q[0];
sx q[0];
rz(8.6046594) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(0.63280025) q[1];
sx q[1];
rz(11.735698) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75738534) q[0];
sx q[0];
rz(-1.877458) q[0];
sx q[0];
rz(-0.39461179) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73859282) q[2];
sx q[2];
rz(-1.2054218) q[2];
sx q[2];
rz(1.5942758) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1066061) q[1];
sx q[1];
rz(-0.88284661) q[1];
sx q[1];
rz(-2.7290542) q[1];
x q[2];
rz(0.72676267) q[3];
sx q[3];
rz(-2.5708377) q[3];
sx q[3];
rz(-1.6970413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43710199) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(3.0214156) q[2];
rz(-1.1581356) q[3];
sx q[3];
rz(-1.3985671) q[3];
sx q[3];
rz(-2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08081089) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(0.91180116) q[0];
rz(-0.78951019) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(-0.3266913) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3347496) q[0];
sx q[0];
rz(-0.93870367) q[0];
sx q[0];
rz(-2.500781) q[0];
rz(0.15906449) q[2];
sx q[2];
rz(-1.1317562) q[2];
sx q[2];
rz(-2.1829) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5527199) q[1];
sx q[1];
rz(-2.1761804) q[1];
sx q[1];
rz(2.6436716) q[1];
x q[2];
rz(-0.88300206) q[3];
sx q[3];
rz(-1.6371173) q[3];
sx q[3];
rz(1.6460713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6313173) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(-1.2871683) q[2];
rz(-0.10989799) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(1.3818285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.592955) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(2.8116995) q[0];
rz(-2.864481) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(2.0842016) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1304504) q[0];
sx q[0];
rz(-1.5913977) q[0];
sx q[0];
rz(-1.9111454) q[0];
rz(-0.34122841) q[2];
sx q[2];
rz(-1.7592906) q[2];
sx q[2];
rz(0.54268062) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.33144618) q[1];
sx q[1];
rz(-0.76008893) q[1];
sx q[1];
rz(1.0815094) q[1];
x q[2];
rz(1.7626761) q[3];
sx q[3];
rz(-1.1550316) q[3];
sx q[3];
rz(-0.92431812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7827591) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(1.3624297) q[2];
rz(-0.6247012) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(-2.3220298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3574922) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(0.5350565) q[0];
rz(2.0013981) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-0.16539703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9388158) q[0];
sx q[0];
rz(-1.6949777) q[0];
sx q[0];
rz(-1.3846272) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.303316) q[2];
sx q[2];
rz(-0.24176134) q[2];
sx q[2];
rz(1.4571112) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0927825) q[1];
sx q[1];
rz(-0.65199344) q[1];
sx q[1];
rz(-0.55744967) q[1];
rz(-pi) q[2];
rz(0.27407077) q[3];
sx q[3];
rz(-1.032853) q[3];
sx q[3];
rz(-0.9243954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.39607221) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(2.9361434) q[2];
rz(-1.127634) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(-1.8849461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(0.35476312) q[0];
rz(-1.154249) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(2.9096471) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1922798) q[0];
sx q[0];
rz(-2.4711907) q[0];
sx q[0];
rz(0.85286661) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7115713) q[2];
sx q[2];
rz(-2.2614334) q[2];
sx q[2];
rz(-1.3123133) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8307454) q[1];
sx q[1];
rz(-1.7880882) q[1];
sx q[1];
rz(-0.22549916) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0911345) q[3];
sx q[3];
rz(-2.7726463) q[3];
sx q[3];
rz(-1.7794533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.140124) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(1.2403963) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(-1.3034472) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0444788) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(-0.98908201) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(-2.1441377) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049012262) q[0];
sx q[0];
rz(-2.2708587) q[0];
sx q[0];
rz(2.9655365) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73156725) q[2];
sx q[2];
rz(-1.288207) q[2];
sx q[2];
rz(-0.76991316) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84038489) q[1];
sx q[1];
rz(-1.9921229) q[1];
sx q[1];
rz(-0.2627443) q[1];
x q[2];
rz(-2.698425) q[3];
sx q[3];
rz(-0.89649761) q[3];
sx q[3];
rz(-0.5815732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.22770195) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(-2.690199) q[2];
rz(0.40870062) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(1.9394978) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8354427) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(2.5174482) q[0];
rz(-1.5165326) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(-2.5278032) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3261624) q[0];
sx q[0];
rz(-0.76061941) q[0];
sx q[0];
rz(-1.8421696) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28366144) q[2];
sx q[2];
rz(-1.5093056) q[2];
sx q[2];
rz(-2.9279857) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2144199) q[1];
sx q[1];
rz(-1.5121636) q[1];
sx q[1];
rz(-0.10417948) q[1];
x q[2];
rz(1.3404487) q[3];
sx q[3];
rz(-0.47603546) q[3];
sx q[3];
rz(-1.7649094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43341407) q[2];
sx q[2];
rz(-0.91579473) q[2];
sx q[2];
rz(1.1748574) q[2];
rz(-2.5332149) q[3];
sx q[3];
rz(-1.7374246) q[3];
sx q[3];
rz(1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90010086) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-2.6532145) q[0];
rz(-1.6237367) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(0.98446313) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0122089) q[0];
sx q[0];
rz(-1.0249656) q[0];
sx q[0];
rz(-1.6145541) q[0];
rz(-2.7266399) q[2];
sx q[2];
rz(-2.3256362) q[2];
sx q[2];
rz(0.50843898) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22673785) q[1];
sx q[1];
rz(-1.2054772) q[1];
sx q[1];
rz(-1.564333) q[1];
x q[2];
rz(0.55926178) q[3];
sx q[3];
rz(-0.54702938) q[3];
sx q[3];
rz(2.3032041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4349334) q[2];
sx q[2];
rz(-1.6985396) q[2];
sx q[2];
rz(-0.53517503) q[2];
rz(1.0501856) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64131367) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(2.4556659) q[0];
rz(-0.39086875) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(-0.92591441) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029322421) q[0];
sx q[0];
rz(-1.7824714) q[0];
sx q[0];
rz(-2.9979343) q[0];
rz(2.4326153) q[2];
sx q[2];
rz(-2.211494) q[2];
sx q[2];
rz(-2.7880653) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8402108) q[1];
sx q[1];
rz(-1.7173319) q[1];
sx q[1];
rz(1.128003) q[1];
rz(-pi) q[2];
rz(2.6191606) q[3];
sx q[3];
rz(-0.2345095) q[3];
sx q[3];
rz(-1.0684551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9459076) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(2.6055028) q[2];
rz(2.7219971) q[3];
sx q[3];
rz(-0.59046888) q[3];
sx q[3];
rz(1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0855899) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(0.39500239) q[0];
rz(-1.5123873) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(1.013247) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7618338) q[0];
sx q[0];
rz(-1.8031617) q[0];
sx q[0];
rz(-2.9985715) q[0];
rz(-pi) q[1];
rz(-0.11864885) q[2];
sx q[2];
rz(-2.7855706) q[2];
sx q[2];
rz(0.3686541) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0777178) q[1];
sx q[1];
rz(-2.7091654) q[1];
sx q[1];
rz(2.9090911) q[1];
rz(-pi) q[2];
rz(-2.4471531) q[3];
sx q[3];
rz(-0.89530066) q[3];
sx q[3];
rz(1.465786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.44832486) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(-2.3821793) q[2];
rz(-1.365472) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7286745) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(2.5683174) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(1.7309932) q[2];
sx q[2];
rz(-0.44993958) q[2];
sx q[2];
rz(-1.1100563) q[2];
rz(-1.6462973) q[3];
sx q[3];
rz(-0.69926881) q[3];
sx q[3];
rz(0.26509501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];