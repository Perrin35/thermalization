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
rz(1.3697019) q[0];
sx q[0];
rz(5.852795) q[0];
sx q[0];
rz(9.2601321) q[0];
rz(1.4511664) q[1];
sx q[1];
rz(-2.7632406) q[1];
sx q[1];
rz(-0.70986706) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71043308) q[0];
sx q[0];
rz(-1.3069777) q[0];
sx q[0];
rz(0.13508787) q[0];
rz(-2.4455382) q[2];
sx q[2];
rz(-2.5230319) q[2];
sx q[2];
rz(1.5355009) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8376417) q[1];
sx q[1];
rz(-1.6198938) q[1];
sx q[1];
rz(1.2173247) q[1];
rz(-2.6534901) q[3];
sx q[3];
rz(-0.26345601) q[3];
sx q[3];
rz(1.7149705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6115438) q[2];
sx q[2];
rz(-0.350746) q[2];
sx q[2];
rz(-1.2190399) q[2];
rz(2.6576095) q[3];
sx q[3];
rz(-0.84570208) q[3];
sx q[3];
rz(-1.774196) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4265863) q[0];
sx q[0];
rz(-0.33892092) q[0];
sx q[0];
rz(1.4434848) q[0];
rz(1.8679856) q[1];
sx q[1];
rz(-2.1111635) q[1];
sx q[1];
rz(2.4748306) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0242321) q[0];
sx q[0];
rz(-1.4015645) q[0];
sx q[0];
rz(0.30099543) q[0];
rz(-pi) q[1];
rz(1.4589492) q[2];
sx q[2];
rz(-1.4750655) q[2];
sx q[2];
rz(-2.7494631) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1135865) q[1];
sx q[1];
rz(-1.1328814) q[1];
sx q[1];
rz(0.44930812) q[1];
x q[2];
rz(-1.9513034) q[3];
sx q[3];
rz(-0.78808053) q[3];
sx q[3];
rz(-2.5684367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.71191177) q[2];
sx q[2];
rz(-1.1744171) q[2];
sx q[2];
rz(-2.5453117) q[2];
rz(1.0239673) q[3];
sx q[3];
rz(-1.7496611) q[3];
sx q[3];
rz(2.021324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4100818) q[0];
sx q[0];
rz(-0.54786587) q[0];
sx q[0];
rz(-1.4672853) q[0];
rz(-0.038854988) q[1];
sx q[1];
rz(-1.7170186) q[1];
sx q[1];
rz(-2.255596) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9666834) q[0];
sx q[0];
rz(-2.0916478) q[0];
sx q[0];
rz(-1.0731927) q[0];
x q[1];
rz(1.6396693) q[2];
sx q[2];
rz(-2.4703272) q[2];
sx q[2];
rz(1.6948989) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3275324) q[1];
sx q[1];
rz(-0.36015727) q[1];
sx q[1];
rz(-3.0292757) q[1];
rz(-2.233611) q[3];
sx q[3];
rz(-3.040979) q[3];
sx q[3];
rz(0.20357547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92750612) q[2];
sx q[2];
rz(-1.4518041) q[2];
sx q[2];
rz(2.5778594) q[2];
rz(-0.8756513) q[3];
sx q[3];
rz(-1.0522269) q[3];
sx q[3];
rz(-2.0503069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57986528) q[0];
sx q[0];
rz(-2.2932597) q[0];
sx q[0];
rz(-2.0735829) q[0];
rz(0.39455286) q[1];
sx q[1];
rz(-2.6119472) q[1];
sx q[1];
rz(-1.4716757) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1009586) q[0];
sx q[0];
rz(-1.9143595) q[0];
sx q[0];
rz(2.6885607) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0288122) q[2];
sx q[2];
rz(-1.28714) q[2];
sx q[2];
rz(0.04087651) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2100239) q[1];
sx q[1];
rz(-1.9940901) q[1];
sx q[1];
rz(-0.84454899) q[1];
rz(-pi) q[2];
rz(1.4626462) q[3];
sx q[3];
rz(-1.9675072) q[3];
sx q[3];
rz(-1.0856398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39845211) q[2];
sx q[2];
rz(-1.8299711) q[2];
sx q[2];
rz(3.0121646) q[2];
rz(0.5433003) q[3];
sx q[3];
rz(-1.0933417) q[3];
sx q[3];
rz(1.3885385) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1035136) q[0];
sx q[0];
rz(-1.0032126) q[0];
sx q[0];
rz(-2.6858618) q[0];
rz(-2.575846) q[1];
sx q[1];
rz(-2.5028298) q[1];
sx q[1];
rz(-2.9260213) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25260293) q[0];
sx q[0];
rz(-0.71160331) q[0];
sx q[0];
rz(1.8322741) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2224765) q[2];
sx q[2];
rz(-1.0813776) q[2];
sx q[2];
rz(0.31351837) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6681155) q[1];
sx q[1];
rz(-1.6398506) q[1];
sx q[1];
rz(1.3051239) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2858538) q[3];
sx q[3];
rz(-1.4654034) q[3];
sx q[3];
rz(1.2440497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4244708) q[2];
sx q[2];
rz(-0.30693808) q[2];
sx q[2];
rz(-1.7090939) q[2];
rz(-1.1411544) q[3];
sx q[3];
rz(-1.2204095) q[3];
sx q[3];
rz(2.6728163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6345217) q[0];
sx q[0];
rz(-2.2880726) q[0];
sx q[0];
rz(-2.5216907) q[0];
rz(-1.9339405) q[1];
sx q[1];
rz(-0.67283216) q[1];
sx q[1];
rz(2.8501453) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7961162) q[0];
sx q[0];
rz(-1.7347851) q[0];
sx q[0];
rz(-1.0046474) q[0];
x q[1];
rz(1.7038149) q[2];
sx q[2];
rz(-1.2499193) q[2];
sx q[2];
rz(1.4747628) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.92293948) q[1];
sx q[1];
rz(-0.99183741) q[1];
sx q[1];
rz(-1.283899) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5813555) q[3];
sx q[3];
rz(-0.91411763) q[3];
sx q[3];
rz(-2.9845723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42346272) q[2];
sx q[2];
rz(-0.66897696) q[2];
sx q[2];
rz(2.1675229) q[2];
rz(-1.2232716) q[3];
sx q[3];
rz(-1.4281102) q[3];
sx q[3];
rz(2.1035002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28737268) q[0];
sx q[0];
rz(-1.1690451) q[0];
sx q[0];
rz(0.30782345) q[0];
rz(2.8616915) q[1];
sx q[1];
rz(-2.8925536) q[1];
sx q[1];
rz(-1.261796) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4331524) q[0];
sx q[0];
rz(-0.13841329) q[0];
sx q[0];
rz(-1.8192151) q[0];
rz(-2.4096411) q[2];
sx q[2];
rz(-1.8667272) q[2];
sx q[2];
rz(-1.0481121) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3615001) q[1];
sx q[1];
rz(-2.8463125) q[1];
sx q[1];
rz(1.5321945) q[1];
x q[2];
rz(-1.7097394) q[3];
sx q[3];
rz(-0.46907237) q[3];
sx q[3];
rz(-1.539849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5032924) q[2];
sx q[2];
rz(-1.1887447) q[2];
sx q[2];
rz(0.80257455) q[2];
rz(-1.59683) q[3];
sx q[3];
rz(-0.39512008) q[3];
sx q[3];
rz(-1.4748658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751727) q[0];
sx q[0];
rz(-1.5479227) q[0];
sx q[0];
rz(-2.5286034) q[0];
rz(2.0892443) q[1];
sx q[1];
rz(-2.3520825) q[1];
sx q[1];
rz(1.8863511) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8156453) q[0];
sx q[0];
rz(-1.1425352) q[0];
sx q[0];
rz(-2.0518579) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1034715) q[2];
sx q[2];
rz(-1.7591068) q[2];
sx q[2];
rz(-1.7323158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71129942) q[1];
sx q[1];
rz(-1.0984527) q[1];
sx q[1];
rz(-0.76763726) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6735821) q[3];
sx q[3];
rz(-2.8492894) q[3];
sx q[3];
rz(1.7287031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.5959979) q[2];
sx q[2];
rz(-0.91384807) q[2];
sx q[2];
rz(0.14860281) q[2];
rz(-0.654486) q[3];
sx q[3];
rz(-1.196967) q[3];
sx q[3];
rz(1.7053846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8504836) q[0];
sx q[0];
rz(-1.608526) q[0];
sx q[0];
rz(-3.1403551) q[0];
rz(1.2507863) q[1];
sx q[1];
rz(-2.2092399) q[1];
sx q[1];
rz(0.95157448) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5244999) q[0];
sx q[0];
rz(-1.6902335) q[0];
sx q[0];
rz(0.16908428) q[0];
rz(-pi) q[1];
rz(-0.96828332) q[2];
sx q[2];
rz(-2.4370425) q[2];
sx q[2];
rz(-0.3144484) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9597541) q[1];
sx q[1];
rz(-1.6153471) q[1];
sx q[1];
rz(-0.049330508) q[1];
rz(-0.22208993) q[3];
sx q[3];
rz(-2.7324893) q[3];
sx q[3];
rz(-0.7717742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9963659) q[2];
sx q[2];
rz(-1.9731015) q[2];
sx q[2];
rz(-2.1053947) q[2];
rz(2.2717617) q[3];
sx q[3];
rz(-1.6706322) q[3];
sx q[3];
rz(-0.74365348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0319801) q[0];
sx q[0];
rz(-0.59015048) q[0];
sx q[0];
rz(-2.0330644) q[0];
rz(-2.176586) q[1];
sx q[1];
rz(-1.3771649) q[1];
sx q[1];
rz(0.56934294) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66273615) q[0];
sx q[0];
rz(-1.3532188) q[0];
sx q[0];
rz(2.6170501) q[0];
x q[1];
rz(-0.37306786) q[2];
sx q[2];
rz(-2.0072492) q[2];
sx q[2];
rz(1.2508357) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7637304) q[1];
sx q[1];
rz(-1.5790577) q[1];
sx q[1];
rz(-2.0935835) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9657992) q[3];
sx q[3];
rz(-2.3758278) q[3];
sx q[3];
rz(2.898223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0137332) q[2];
sx q[2];
rz(-0.49488417) q[2];
sx q[2];
rz(-0.046517046) q[2];
rz(-2.8401996) q[3];
sx q[3];
rz(-1.2208341) q[3];
sx q[3];
rz(2.9197599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(2.2409869) q[0];
sx q[0];
rz(-1.2099246) q[0];
sx q[0];
rz(-1.8151617) q[0];
rz(1.2251414) q[1];
sx q[1];
rz(-0.50347181) q[1];
sx q[1];
rz(1.0540963) q[1];
rz(2.8711748) q[2];
sx q[2];
rz(-2.5106988) q[2];
sx q[2];
rz(-0.0028263447) q[2];
rz(-0.75792652) q[3];
sx q[3];
rz(-1.0316385) q[3];
sx q[3];
rz(1.3760174) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
