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
rz(-1.2404233) q[0];
sx q[0];
rz(-1.197553) q[0];
sx q[0];
rz(2.9252606) q[0];
rz(4.0989838) q[1];
sx q[1];
rz(5.6070072) q[1];
sx q[1];
rz(11.054872) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0441598) q[0];
sx q[0];
rz(-0.55632797) q[0];
sx q[0];
rz(-3.1233643) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47343238) q[2];
sx q[2];
rz(-1.2780927) q[2];
sx q[2];
rz(2.6930489) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2701157) q[1];
sx q[1];
rz(-1.9028712) q[1];
sx q[1];
rz(0.52712743) q[1];
rz(-pi) q[2];
rz(1.1153631) q[3];
sx q[3];
rz(-0.75616992) q[3];
sx q[3];
rz(2.2952473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.034885255) q[2];
sx q[2];
rz(-0.35197508) q[2];
sx q[2];
rz(1.0484288) q[2];
rz(0.18167051) q[3];
sx q[3];
rz(-2.1766267) q[3];
sx q[3];
rz(2.488193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.492391) q[0];
sx q[0];
rz(-0.94434706) q[0];
sx q[0];
rz(2.6948068) q[0];
rz(2.2611639) q[1];
sx q[1];
rz(-1.7767521) q[1];
sx q[1];
rz(0.78278881) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.370011) q[0];
sx q[0];
rz(-0.25213045) q[0];
sx q[0];
rz(1.3229516) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1433463) q[2];
sx q[2];
rz(-0.99787092) q[2];
sx q[2];
rz(-1.7597511) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.782541) q[1];
sx q[1];
rz(-2.5757669) q[1];
sx q[1];
rz(0.58360175) q[1];
rz(0.47329013) q[3];
sx q[3];
rz(-0.9447228) q[3];
sx q[3];
rz(-1.9260709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.11144) q[2];
sx q[2];
rz(-1.4806662) q[2];
sx q[2];
rz(2.1622369) q[2];
rz(2.7566946) q[3];
sx q[3];
rz(-1.9210457) q[3];
sx q[3];
rz(2.9773007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48591831) q[0];
sx q[0];
rz(-3.0874708) q[0];
sx q[0];
rz(0.78980494) q[0];
rz(0.18665953) q[1];
sx q[1];
rz(-1.7402382) q[1];
sx q[1];
rz(-0.99367118) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45897608) q[0];
sx q[0];
rz(-2.1135215) q[0];
sx q[0];
rz(-1.2786464) q[0];
rz(-pi) q[1];
rz(0.61305253) q[2];
sx q[2];
rz(-1.7568577) q[2];
sx q[2];
rz(1.8579872) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.51998617) q[1];
sx q[1];
rz(-2.4944759) q[1];
sx q[1];
rz(-2.3220329) q[1];
rz(1.8422115) q[3];
sx q[3];
rz(-0.86835734) q[3];
sx q[3];
rz(-2.0347119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32187244) q[2];
sx q[2];
rz(-1.1178144) q[2];
sx q[2];
rz(-0.37128386) q[2];
rz(-2.7927981) q[3];
sx q[3];
rz(-1.0943509) q[3];
sx q[3];
rz(0.5955407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6868941) q[0];
sx q[0];
rz(-1.0201539) q[0];
sx q[0];
rz(1.4042847) q[0];
rz(0.67717254) q[1];
sx q[1];
rz(-1.155747) q[1];
sx q[1];
rz(-1.6489395) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9335404) q[0];
sx q[0];
rz(-1.786199) q[0];
sx q[0];
rz(-0.24788863) q[0];
x q[1];
rz(-2.1588232) q[2];
sx q[2];
rz(-0.50543907) q[2];
sx q[2];
rz(2.7014521) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4863534) q[1];
sx q[1];
rz(-0.35105536) q[1];
sx q[1];
rz(0.69610657) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2954313) q[3];
sx q[3];
rz(-1.3337047) q[3];
sx q[3];
rz(-0.8704291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6877785) q[2];
sx q[2];
rz(-0.9684338) q[2];
sx q[2];
rz(0.34823927) q[2];
rz(1.4549152) q[3];
sx q[3];
rz(-1.429052) q[3];
sx q[3];
rz(-1.7961563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1763879) q[0];
sx q[0];
rz(-2.5768953) q[0];
sx q[0];
rz(-2.2494466) q[0];
rz(-2.6773894) q[1];
sx q[1];
rz(-1.2449539) q[1];
sx q[1];
rz(1.8468599) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2703122) q[0];
sx q[0];
rz(-1.082232) q[0];
sx q[0];
rz(-1.0315007) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6038997) q[2];
sx q[2];
rz(-1.6811996) q[2];
sx q[2];
rz(2.3619719) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6695822) q[1];
sx q[1];
rz(-1.246844) q[1];
sx q[1];
rz(-2.0355909) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61035778) q[3];
sx q[3];
rz(-0.98516432) q[3];
sx q[3];
rz(2.741339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3225473) q[2];
sx q[2];
rz(-0.61083856) q[2];
sx q[2];
rz(-0.56582212) q[2];
rz(0.081341751) q[3];
sx q[3];
rz(-2.1606725) q[3];
sx q[3];
rz(2.5210023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.464798) q[0];
sx q[0];
rz(-3.0749574) q[0];
sx q[0];
rz(-1.5555405) q[0];
rz(2.0809035) q[1];
sx q[1];
rz(-1.5763177) q[1];
sx q[1];
rz(-0.63180822) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2151129) q[0];
sx q[0];
rz(-0.2479015) q[0];
sx q[0];
rz(0.32899022) q[0];
rz(-pi) q[1];
rz(-1.00497) q[2];
sx q[2];
rz(-0.7938876) q[2];
sx q[2];
rz(-0.53295202) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.28367701) q[1];
sx q[1];
rz(-1.9092968) q[1];
sx q[1];
rz(-1.0033016) q[1];
x q[2];
rz(-1.7177714) q[3];
sx q[3];
rz(-0.89255652) q[3];
sx q[3];
rz(3.105046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1598728) q[2];
sx q[2];
rz(-2.0231415) q[2];
sx q[2];
rz(2.6427606) q[2];
rz(1.3129129) q[3];
sx q[3];
rz(-2.4098318) q[3];
sx q[3];
rz(1.7074728) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6071534) q[0];
sx q[0];
rz(-2.7767015) q[0];
sx q[0];
rz(0.69951192) q[0];
rz(2.6761159) q[1];
sx q[1];
rz(-0.87124467) q[1];
sx q[1];
rz(-2.2043998) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4995183) q[0];
sx q[0];
rz(-2.8694911) q[0];
sx q[0];
rz(0.67263747) q[0];
rz(-pi) q[1];
rz(2.7363051) q[2];
sx q[2];
rz(-2.8817085) q[2];
sx q[2];
rz(-1.056162) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73884642) q[1];
sx q[1];
rz(-2.7838696) q[1];
sx q[1];
rz(0.92836942) q[1];
x q[2];
rz(-2.3714957) q[3];
sx q[3];
rz(-0.40989629) q[3];
sx q[3];
rz(2.485825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.84404868) q[2];
sx q[2];
rz(-1.301845) q[2];
sx q[2];
rz(0.37332264) q[2];
rz(-1.9518055) q[3];
sx q[3];
rz(-2.6326284) q[3];
sx q[3];
rz(2.6313307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.74476403) q[0];
sx q[0];
rz(-1.0823534) q[0];
sx q[0];
rz(0.74791351) q[0];
rz(0.76639908) q[1];
sx q[1];
rz(-0.26873573) q[1];
sx q[1];
rz(3.1386197) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.317694) q[0];
sx q[0];
rz(-2.4269773) q[0];
sx q[0];
rz(1.8665642) q[0];
rz(-pi) q[1];
rz(3.1220436) q[2];
sx q[2];
rz(-1.4998933) q[2];
sx q[2];
rz(-0.30724684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.06719477) q[1];
sx q[1];
rz(-1.7247827) q[1];
sx q[1];
rz(-0.043937307) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3467196) q[3];
sx q[3];
rz(-2.8668781) q[3];
sx q[3];
rz(2.2661254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.030674) q[2];
sx q[2];
rz(-1.3838394) q[2];
sx q[2];
rz(1.0670916) q[2];
rz(0.083960697) q[3];
sx q[3];
rz(-2.6559918) q[3];
sx q[3];
rz(-2.3626204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.344051) q[0];
sx q[0];
rz(-2.2028956) q[0];
sx q[0];
rz(3.0352266) q[0];
rz(-2.1616409) q[1];
sx q[1];
rz(-1.4915024) q[1];
sx q[1];
rz(2.3756557) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.879012) q[0];
sx q[0];
rz(-2.0011138) q[0];
sx q[0];
rz(0.58027123) q[0];
rz(-pi) q[1];
rz(-1.4808876) q[2];
sx q[2];
rz(-1.3718296) q[2];
sx q[2];
rz(-2.112191) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3801743) q[1];
sx q[1];
rz(-1.5360218) q[1];
sx q[1];
rz(-2.4654287) q[1];
rz(0.25772734) q[3];
sx q[3];
rz(-0.014583909) q[3];
sx q[3];
rz(-0.38250438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.47436675) q[2];
sx q[2];
rz(-0.52046481) q[2];
sx q[2];
rz(-0.70029798) q[2];
rz(2.4750366) q[3];
sx q[3];
rz(-1.8066112) q[3];
sx q[3];
rz(-1.0880281) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661082) q[0];
sx q[0];
rz(-1.0270783) q[0];
sx q[0];
rz(-1.745537) q[0];
rz(-1.4736157) q[1];
sx q[1];
rz(-1.8831848) q[1];
sx q[1];
rz(-0.6667164) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5294801) q[0];
sx q[0];
rz(-0.94769883) q[0];
sx q[0];
rz(0.77147958) q[0];
rz(-pi) q[1];
rz(0.36138968) q[2];
sx q[2];
rz(-2.4494684) q[2];
sx q[2];
rz(-3.1049984) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.085070193) q[1];
sx q[1];
rz(-1.4398972) q[1];
sx q[1];
rz(2.8225949) q[1];
x q[2];
rz(2.5693043) q[3];
sx q[3];
rz(-1.0977355) q[3];
sx q[3];
rz(2.8561887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.066976808) q[2];
sx q[2];
rz(-2.4846027) q[2];
sx q[2];
rz(-2.2508049) q[2];
rz(2.7095419) q[3];
sx q[3];
rz(-1.0703577) q[3];
sx q[3];
rz(0.74705684) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73410949) q[0];
sx q[0];
rz(-1.4984087) q[0];
sx q[0];
rz(1.8468504) q[0];
rz(1.2474077) q[1];
sx q[1];
rz(-0.71949646) q[1];
sx q[1];
rz(1.234642) q[1];
rz(0.16130372) q[2];
sx q[2];
rz(-1.9111173) q[2];
sx q[2];
rz(2.5572122) q[2];
rz(-2.6525146) q[3];
sx q[3];
rz(-1.5413956) q[3];
sx q[3];
rz(1.2914381) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
