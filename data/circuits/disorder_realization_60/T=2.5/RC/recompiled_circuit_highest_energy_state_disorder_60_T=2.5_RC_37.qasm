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
rz(-2.5118339) q[0];
sx q[0];
rz(-0.30183733) q[0];
sx q[0];
rz(-1.8344185) q[0];
rz(2.6131926) q[1];
sx q[1];
rz(-2.3101248) q[1];
sx q[1];
rz(2.6822579) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0889657) q[0];
sx q[0];
rz(-1.8267434) q[0];
sx q[0];
rz(3.0908282) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9436753) q[2];
sx q[2];
rz(-1.8423586) q[2];
sx q[2];
rz(2.7261734) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6428549) q[1];
sx q[1];
rz(-0.6371405) q[1];
sx q[1];
rz(-2.8481507) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.022597952) q[3];
sx q[3];
rz(-2.239478) q[3];
sx q[3];
rz(-0.99136664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.022973148) q[2];
sx q[2];
rz(-1.8034673) q[2];
sx q[2];
rz(-2.5845134) q[2];
rz(-2.6856375) q[3];
sx q[3];
rz(-0.7619226) q[3];
sx q[3];
rz(2.5383811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93422455) q[0];
sx q[0];
rz(-0.71894431) q[0];
sx q[0];
rz(-1.7508605) q[0];
rz(-1.8234183) q[1];
sx q[1];
rz(-1.7134106) q[1];
sx q[1];
rz(0.21600977) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0586313) q[0];
sx q[0];
rz(-2.4602232) q[0];
sx q[0];
rz(1.0555223) q[0];
rz(-1.4513218) q[2];
sx q[2];
rz(-1.9873957) q[2];
sx q[2];
rz(-2.4753776) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5358754) q[1];
sx q[1];
rz(-0.89369666) q[1];
sx q[1];
rz(0.89800055) q[1];
rz(-pi) q[2];
rz(-2.1116637) q[3];
sx q[3];
rz(-1.7591624) q[3];
sx q[3];
rz(-2.4595367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7184489) q[2];
sx q[2];
rz(-2.7687912) q[2];
sx q[2];
rz(-0.049840363) q[2];
rz(-0.54346624) q[3];
sx q[3];
rz(-1.1246357) q[3];
sx q[3];
rz(-2.0487093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24930382) q[0];
sx q[0];
rz(-2.0396621) q[0];
sx q[0];
rz(-2.1732543) q[0];
rz(-3.1210506) q[1];
sx q[1];
rz(-1.3803218) q[1];
sx q[1];
rz(1.4302018) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1002625) q[0];
sx q[0];
rz(-1.593129) q[0];
sx q[0];
rz(-2.5002527) q[0];
rz(2.5920141) q[2];
sx q[2];
rz(-2.5039154) q[2];
sx q[2];
rz(-2.9108832) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80950981) q[1];
sx q[1];
rz(-0.44594279) q[1];
sx q[1];
rz(0.46582241) q[1];
rz(-pi) q[2];
rz(-0.3768206) q[3];
sx q[3];
rz(-1.539789) q[3];
sx q[3];
rz(-0.12801192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4358383) q[2];
sx q[2];
rz(-2.4987554) q[2];
sx q[2];
rz(-1.7819116) q[2];
rz(-2.6333574) q[3];
sx q[3];
rz(-1.5225007) q[3];
sx q[3];
rz(-1.9967509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.085676) q[0];
sx q[0];
rz(-2.8450232) q[0];
sx q[0];
rz(-2.1348409) q[0];
rz(-2.309917) q[1];
sx q[1];
rz(-0.74200231) q[1];
sx q[1];
rz(-2.0337909) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.09064085) q[0];
sx q[0];
rz(-1.9371175) q[0];
sx q[0];
rz(-0.90496814) q[0];
x q[1];
rz(2.7659589) q[2];
sx q[2];
rz(-1.8058597) q[2];
sx q[2];
rz(-1.3974853) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9398168) q[1];
sx q[1];
rz(-0.30758938) q[1];
sx q[1];
rz(1.5053476) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8516225) q[3];
sx q[3];
rz(-1.1410332) q[3];
sx q[3];
rz(-2.6971779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.45336777) q[2];
sx q[2];
rz(-1.0658762) q[2];
sx q[2];
rz(-0.31309703) q[2];
rz(-2.744216) q[3];
sx q[3];
rz(-1.671096) q[3];
sx q[3];
rz(0.1853005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.6322286) q[0];
sx q[0];
rz(-2.2780184) q[0];
sx q[0];
rz(-2.2160227) q[0];
rz(2.6717692) q[1];
sx q[1];
rz(-1.8269822) q[1];
sx q[1];
rz(2.125461) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39611577) q[0];
sx q[0];
rz(-1.9821616) q[0];
sx q[0];
rz(-1.0883254) q[0];
rz(-pi) q[1];
rz(2.8723628) q[2];
sx q[2];
rz(-1.2310264) q[2];
sx q[2];
rz(1.4087848) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4279856) q[1];
sx q[1];
rz(-1.3283037) q[1];
sx q[1];
rz(2.9338783) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1096439) q[3];
sx q[3];
rz(-1.2476139) q[3];
sx q[3];
rz(1.1695216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4887345) q[2];
sx q[2];
rz(-2.8159339) q[2];
sx q[2];
rz(-0.48750901) q[2];
rz(2.2440535) q[3];
sx q[3];
rz(-1.2582015) q[3];
sx q[3];
rz(-1.0645617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2010736) q[0];
sx q[0];
rz(-0.93795332) q[0];
sx q[0];
rz(2.4679389) q[0];
rz(-1.9081839) q[1];
sx q[1];
rz(-1.0181095) q[1];
sx q[1];
rz(2.3755551) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4318643) q[0];
sx q[0];
rz(-2.0805295) q[0];
sx q[0];
rz(-2.8317189) q[0];
rz(-pi) q[1];
rz(1.4991708) q[2];
sx q[2];
rz(-1.7461458) q[2];
sx q[2];
rz(-1.2855114) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9711069) q[1];
sx q[1];
rz(-2.0506713) q[1];
sx q[1];
rz(-0.1636774) q[1];
rz(-pi) q[2];
rz(-2.8358342) q[3];
sx q[3];
rz(-1.4396724) q[3];
sx q[3];
rz(2.6763889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.96887863) q[2];
sx q[2];
rz(-1.6438899) q[2];
sx q[2];
rz(2.3806351) q[2];
rz(-2.739665) q[3];
sx q[3];
rz(-2.8765078) q[3];
sx q[3];
rz(-0.5886122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6589979) q[0];
sx q[0];
rz(-2.3880279) q[0];
sx q[0];
rz(1.6931417) q[0];
rz(-2.6407369) q[1];
sx q[1];
rz(-1.294699) q[1];
sx q[1];
rz(-2.3282611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017352176) q[0];
sx q[0];
rz(-0.4158786) q[0];
sx q[0];
rz(-2.7680567) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0221239) q[2];
sx q[2];
rz(-1.0257693) q[2];
sx q[2];
rz(-1.8785005) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.2351027) q[1];
sx q[1];
rz(-1.1336599) q[1];
sx q[1];
rz(1.7227931) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1239979) q[3];
sx q[3];
rz(-0.13544336) q[3];
sx q[3];
rz(2.9970084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8616051) q[2];
sx q[2];
rz(-0.6051175) q[2];
sx q[2];
rz(-2.81847) q[2];
rz(-1.4019639) q[3];
sx q[3];
rz(-1.8596545) q[3];
sx q[3];
rz(-1.3824979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.062227) q[0];
sx q[0];
rz(-2.0923738) q[0];
sx q[0];
rz(2.3754689) q[0];
rz(2.3221723) q[1];
sx q[1];
rz(-1.0304281) q[1];
sx q[1];
rz(1.5310418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9510244) q[0];
sx q[0];
rz(-1.5026717) q[0];
sx q[0];
rz(-0.22916746) q[0];
x q[1];
rz(0.28240164) q[2];
sx q[2];
rz(-0.67495433) q[2];
sx q[2];
rz(-0.16004983) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4447915) q[1];
sx q[1];
rz(-2.2292622) q[1];
sx q[1];
rz(-2.3347003) q[1];
rz(-pi) q[2];
rz(0.29115486) q[3];
sx q[3];
rz(-1.8829573) q[3];
sx q[3];
rz(-0.96033421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3766342) q[2];
sx q[2];
rz(-2.9464293) q[2];
sx q[2];
rz(0.39806077) q[2];
rz(0.6063439) q[3];
sx q[3];
rz(-1.5667934) q[3];
sx q[3];
rz(-2.2280367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9048318) q[0];
sx q[0];
rz(-0.30192152) q[0];
sx q[0];
rz(2.6171369) q[0];
rz(0.66871387) q[1];
sx q[1];
rz(-1.5293744) q[1];
sx q[1];
rz(1.3800157) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83630122) q[0];
sx q[0];
rz(-0.99812767) q[0];
sx q[0];
rz(-0.051866626) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0247714) q[2];
sx q[2];
rz(-2.5971892) q[2];
sx q[2];
rz(2.4449206) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.42889574) q[1];
sx q[1];
rz(-1.98723) q[1];
sx q[1];
rz(0.45559366) q[1];
rz(-pi) q[2];
rz(1.947375) q[3];
sx q[3];
rz(-2.0717952) q[3];
sx q[3];
rz(-0.11160103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1606719) q[2];
sx q[2];
rz(-2.394684) q[2];
sx q[2];
rz(0.46372947) q[2];
rz(1.7175698) q[3];
sx q[3];
rz(-1.0878891) q[3];
sx q[3];
rz(0.033528479) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5793107) q[0];
sx q[0];
rz(-0.71826851) q[0];
sx q[0];
rz(2.723519) q[0];
rz(2.383291) q[1];
sx q[1];
rz(-2.1577991) q[1];
sx q[1];
rz(-1.2850579) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4649974) q[0];
sx q[0];
rz(-0.9228068) q[0];
sx q[0];
rz(1.4761488) q[0];
rz(-pi) q[1];
rz(2.7964716) q[2];
sx q[2];
rz(-2.1639544) q[2];
sx q[2];
rz(-2.5760004) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1753208) q[1];
sx q[1];
rz(-2.113225) q[1];
sx q[1];
rz(1.7955417) q[1];
rz(-pi) q[2];
rz(-1.9081014) q[3];
sx q[3];
rz(-1.3618338) q[3];
sx q[3];
rz(-1.8551872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.16984223) q[2];
sx q[2];
rz(-2.1368133) q[2];
sx q[2];
rz(0.24610914) q[2];
rz(-1.2717815) q[3];
sx q[3];
rz(-1.5747986) q[3];
sx q[3];
rz(-1.3625712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1590189) q[0];
sx q[0];
rz(-1.6652501) q[0];
sx q[0];
rz(-0.2022947) q[0];
rz(-2.5166439) q[1];
sx q[1];
rz(-0.85660558) q[1];
sx q[1];
rz(0.4013335) q[1];
rz(-1.3622147) q[2];
sx q[2];
rz(-1.3857532) q[2];
sx q[2];
rz(-2.9278853) q[2];
rz(-0.20062867) q[3];
sx q[3];
rz(-1.8798141) q[3];
sx q[3];
rz(0.78165913) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
