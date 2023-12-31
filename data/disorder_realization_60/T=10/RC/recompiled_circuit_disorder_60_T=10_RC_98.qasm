OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1922167) q[0];
sx q[0];
rz(-2.0944216) q[0];
sx q[0];
rz(3.0728683) q[0];
rz(1.7460495) q[1];
sx q[1];
rz(4.6739251) q[1];
sx q[1];
rz(8.2164017) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8954593) q[0];
sx q[0];
rz(-1.6172505) q[0];
sx q[0];
rz(1.6669271) q[0];
rz(-3.0460998) q[2];
sx q[2];
rz(-3.0152233) q[2];
sx q[2];
rz(-0.40590826) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5282643) q[1];
sx q[1];
rz(-1.4600888) q[1];
sx q[1];
rz(-0.61943357) q[1];
rz(-0.73322202) q[3];
sx q[3];
rz(-2.4823722) q[3];
sx q[3];
rz(-2.1472907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8157114) q[2];
sx q[2];
rz(-1.7408966) q[2];
sx q[2];
rz(1.4665843) q[2];
rz(-0.69774929) q[3];
sx q[3];
rz(-2.0402699) q[3];
sx q[3];
rz(2.3944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-3.0242457) q[0];
sx q[0];
rz(-1.1276561) q[0];
sx q[0];
rz(1.1741937) q[0];
rz(-0.17114561) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(-0.29719621) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78404616) q[0];
sx q[0];
rz(-1.6539126) q[0];
sx q[0];
rz(-3.0931285) q[0];
x q[1];
rz(-2.6058795) q[2];
sx q[2];
rz(-2.3539054) q[2];
sx q[2];
rz(0.96131575) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.17644037) q[1];
sx q[1];
rz(-2.0808176) q[1];
sx q[1];
rz(0.39217197) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7278565) q[3];
sx q[3];
rz(-2.3351151) q[3];
sx q[3];
rz(-1.8959351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.16195665) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(0.61398181) q[2];
rz(-0.87614122) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(0.63703018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0959594) q[0];
sx q[0];
rz(-1.8563844) q[0];
sx q[0];
rz(0.30763787) q[0];
rz(0.74854198) q[1];
sx q[1];
rz(-0.33154878) q[1];
sx q[1];
rz(2.3017853) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9687942) q[0];
sx q[0];
rz(-1.8546687) q[0];
sx q[0];
rz(2.7127405) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8790881) q[2];
sx q[2];
rz(-2.7911107) q[2];
sx q[2];
rz(-2.1413435) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0963904) q[1];
sx q[1];
rz(-0.55711105) q[1];
sx q[1];
rz(-0.60369173) q[1];
x q[2];
rz(1.3358467) q[3];
sx q[3];
rz(-2.4953105) q[3];
sx q[3];
rz(0.39138734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86747375) q[2];
sx q[2];
rz(-1.7904736) q[2];
sx q[2];
rz(-1.8236558) q[2];
rz(-1.2157724) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(-1.6320451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3777305) q[0];
sx q[0];
rz(-1.7611935) q[0];
sx q[0];
rz(0.44556251) q[0];
rz(0.62082779) q[1];
sx q[1];
rz(-1.2460243) q[1];
sx q[1];
rz(-0.96558085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27698101) q[0];
sx q[0];
rz(-2.3267713) q[0];
sx q[0];
rz(0.38397249) q[0];
rz(3.1104286) q[2];
sx q[2];
rz(-2.2989797) q[2];
sx q[2];
rz(0.63809168) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.671215) q[1];
sx q[1];
rz(-2.885985) q[1];
sx q[1];
rz(1.6971991) q[1];
x q[2];
rz(-1.6846893) q[3];
sx q[3];
rz(-1.6674862) q[3];
sx q[3];
rz(0.25988042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3199557) q[2];
sx q[2];
rz(-2.379202) q[2];
sx q[2];
rz(0.92932534) q[2];
rz(-0.6435414) q[3];
sx q[3];
rz(-2.1054335) q[3];
sx q[3];
rz(2.1381366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4716361) q[0];
sx q[0];
rz(-1.5595373) q[0];
sx q[0];
rz(-1.8575645) q[0];
rz(-0.28981003) q[1];
sx q[1];
rz(-2.402014) q[1];
sx q[1];
rz(1.0481542) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50131932) q[0];
sx q[0];
rz(-0.87448705) q[0];
sx q[0];
rz(-2.2107844) q[0];
rz(-pi) q[1];
rz(-2.5571312) q[2];
sx q[2];
rz(-2.2197154) q[2];
sx q[2];
rz(1.4635758) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.84761274) q[1];
sx q[1];
rz(-1.8925397) q[1];
sx q[1];
rz(-1.4732248) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9976451) q[3];
sx q[3];
rz(-0.76612681) q[3];
sx q[3];
rz(-1.0539953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8706878) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(-2.2407545) q[2];
rz(-2.0488996) q[3];
sx q[3];
rz(-1.0034424) q[3];
sx q[3];
rz(-1.2341011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1822405) q[0];
sx q[0];
rz(-1.933796) q[0];
sx q[0];
rz(-2.545488) q[0];
rz(-1.6456564) q[1];
sx q[1];
rz(-0.89769617) q[1];
sx q[1];
rz(1.2449107) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7398864) q[0];
sx q[0];
rz(-0.6131999) q[0];
sx q[0];
rz(-0.99341157) q[0];
rz(-pi) q[1];
rz(-2.5858324) q[2];
sx q[2];
rz(-0.55609497) q[2];
sx q[2];
rz(-0.13242002) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.54904304) q[1];
sx q[1];
rz(-0.92714308) q[1];
sx q[1];
rz(-2.2085269) q[1];
rz(0.75002589) q[3];
sx q[3];
rz(-0.3332899) q[3];
sx q[3];
rz(2.1961574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9101377) q[2];
sx q[2];
rz(-2.4427876) q[2];
sx q[2];
rz(0.22496741) q[2];
rz(-3.0531626) q[3];
sx q[3];
rz(-1.6902573) q[3];
sx q[3];
rz(0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6769619) q[0];
sx q[0];
rz(-1.9043652) q[0];
sx q[0];
rz(0.27994573) q[0];
rz(-1.6784558) q[1];
sx q[1];
rz(-1.8755553) q[1];
sx q[1];
rz(-2.8889012) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8425927) q[0];
sx q[0];
rz(-0.29688603) q[0];
sx q[0];
rz(1.7345558) q[0];
rz(-pi) q[1];
rz(-1.5451317) q[2];
sx q[2];
rz(-2.5070094) q[2];
sx q[2];
rz(2.6013825) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.60755) q[1];
sx q[1];
rz(-2.0760963) q[1];
sx q[1];
rz(-2.3984548) q[1];
rz(-pi) q[2];
x q[2];
rz(2.383109) q[3];
sx q[3];
rz(-1.0060203) q[3];
sx q[3];
rz(1.5837216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.32020405) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(-1.6097216) q[2];
rz(-1.1931217) q[3];
sx q[3];
rz(-0.90819287) q[3];
sx q[3];
rz(-2.0549205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10483345) q[0];
sx q[0];
rz(-1.3830673) q[0];
sx q[0];
rz(-1.4861134) q[0];
rz(0.2688109) q[1];
sx q[1];
rz(-2.0116282) q[1];
sx q[1];
rz(2.862646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4501707) q[0];
sx q[0];
rz(-2.2922463) q[0];
sx q[0];
rz(-0.88766092) q[0];
rz(-pi) q[1];
rz(2.9697044) q[2];
sx q[2];
rz(-1.0749146) q[2];
sx q[2];
rz(1.9672729) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3153509) q[1];
sx q[1];
rz(-0.8289753) q[1];
sx q[1];
rz(-0.51893236) q[1];
rz(-0.14848498) q[3];
sx q[3];
rz(-1.0035702) q[3];
sx q[3];
rz(0.70311577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.802861) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(1.5926682) q[2];
rz(2.0907949) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(-1.1999493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46362296) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(-1.8883702) q[0];
rz(-2.5190917) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(-1.1463096) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63908731) q[0];
sx q[0];
rz(-0.82477942) q[0];
sx q[0];
rz(-1.6665002) q[0];
rz(-pi) q[1];
rz(-0.65593221) q[2];
sx q[2];
rz(-2.1428875) q[2];
sx q[2];
rz(-0.42665542) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1100626) q[1];
sx q[1];
rz(-1.4974125) q[1];
sx q[1];
rz(0.88295464) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9911489) q[3];
sx q[3];
rz(-2.162809) q[3];
sx q[3];
rz(2.1136485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.15381947) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(1.0158319) q[2];
rz(-0.9097957) q[3];
sx q[3];
rz(-2.4715021) q[3];
sx q[3];
rz(2.773496) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1290454) q[0];
sx q[0];
rz(-2.5279901) q[0];
sx q[0];
rz(3.1273499) q[0];
rz(0.8447389) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(1.7600118) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2505014) q[0];
sx q[0];
rz(-1.7133461) q[0];
sx q[0];
rz(3.13091) q[0];
rz(-2.1496885) q[2];
sx q[2];
rz(-1.6245981) q[2];
sx q[2];
rz(2.3827041) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4631008) q[1];
sx q[1];
rz(-1.7144202) q[1];
sx q[1];
rz(-2.1302057) q[1];
rz(2.6081309) q[3];
sx q[3];
rz(-1.3812314) q[3];
sx q[3];
rz(1.5272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0649197) q[2];
sx q[2];
rz(-1.2023456) q[2];
sx q[2];
rz(-0.98999611) q[2];
rz(-2.2475217) q[3];
sx q[3];
rz(-0.48864135) q[3];
sx q[3];
rz(-2.9161684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0902696) q[0];
sx q[0];
rz(-2.4933503) q[0];
sx q[0];
rz(-1.1052263) q[0];
rz(1.3399711) q[1];
sx q[1];
rz(-0.62146386) q[1];
sx q[1];
rz(0.38846831) q[1];
rz(0.71628911) q[2];
sx q[2];
rz(-1.204797) q[2];
sx q[2];
rz(2.7439678) q[2];
rz(-2.7064825) q[3];
sx q[3];
rz(-2.1463487) q[3];
sx q[3];
rz(-1.6456732) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
