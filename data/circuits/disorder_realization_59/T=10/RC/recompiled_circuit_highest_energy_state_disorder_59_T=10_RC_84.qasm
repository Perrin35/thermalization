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
rz(-1.8772797) q[0];
sx q[0];
rz(-0.32353434) q[0];
sx q[0];
rz(2.6927595) q[0];
rz(-1.1558865) q[1];
sx q[1];
rz(-1.7855676) q[1];
sx q[1];
rz(-1.1745656) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4455152) q[0];
sx q[0];
rz(-0.68176523) q[0];
sx q[0];
rz(0.56607874) q[0];
x q[1];
rz(-0.20148142) q[2];
sx q[2];
rz(-0.32989943) q[2];
sx q[2];
rz(-0.13452521) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1303181) q[1];
sx q[1];
rz(-2.5198333) q[1];
sx q[1];
rz(1.51062) q[1];
rz(-pi) q[2];
rz(-2.6277718) q[3];
sx q[3];
rz(-1.658545) q[3];
sx q[3];
rz(0.82753554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.590098) q[2];
sx q[2];
rz(-1.8415201) q[2];
sx q[2];
rz(-2.2287915) q[2];
rz(1.270795) q[3];
sx q[3];
rz(-2.1015034) q[3];
sx q[3];
rz(-0.84996581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0047282334) q[0];
sx q[0];
rz(-0.57924634) q[0];
sx q[0];
rz(-0.42214033) q[0];
rz(2.7719851) q[1];
sx q[1];
rz(-1.0969176) q[1];
sx q[1];
rz(0.94013989) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62806118) q[0];
sx q[0];
rz(-1.9002751) q[0];
sx q[0];
rz(-0.0026436289) q[0];
rz(-pi) q[1];
rz(-1.7707509) q[2];
sx q[2];
rz(-2.5995289) q[2];
sx q[2];
rz(2.0393537) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.247929) q[1];
sx q[1];
rz(-1.0808696) q[1];
sx q[1];
rz(-2.41416) q[1];
x q[2];
rz(0.27249725) q[3];
sx q[3];
rz(-2.2292308) q[3];
sx q[3];
rz(0.30425018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1993572) q[2];
sx q[2];
rz(-0.67023674) q[2];
sx q[2];
rz(-2.2352236) q[2];
rz(2.1622315) q[3];
sx q[3];
rz(-2.7370079) q[3];
sx q[3];
rz(-1.8766859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1449428) q[0];
sx q[0];
rz(-0.066983797) q[0];
sx q[0];
rz(-1.2445194) q[0];
rz(2.7527346) q[1];
sx q[1];
rz(-1.8617947) q[1];
sx q[1];
rz(-0.32274524) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53777107) q[0];
sx q[0];
rz(-2.4706998) q[0];
sx q[0];
rz(-0.18108271) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1006196) q[2];
sx q[2];
rz(-2.5289446) q[2];
sx q[2];
rz(-0.27951159) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1908993) q[1];
sx q[1];
rz(-2.3436693) q[1];
sx q[1];
rz(1.8969593) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6204319) q[3];
sx q[3];
rz(-1.921504) q[3];
sx q[3];
rz(1.8010822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2286223) q[2];
sx q[2];
rz(-1.7614438) q[2];
sx q[2];
rz(-1.7886394) q[2];
rz(0.27255034) q[3];
sx q[3];
rz(-1.291178) q[3];
sx q[3];
rz(1.6779617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6005741) q[0];
sx q[0];
rz(-0.32222846) q[0];
sx q[0];
rz(1.898265) q[0];
rz(-0.86878949) q[1];
sx q[1];
rz(-0.96005762) q[1];
sx q[1];
rz(-1.0702081) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6400498) q[0];
sx q[0];
rz(-2.2794855) q[0];
sx q[0];
rz(-0.88592822) q[0];
x q[1];
rz(-0.44311503) q[2];
sx q[2];
rz(-2.296148) q[2];
sx q[2];
rz(1.0034753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8319885) q[1];
sx q[1];
rz(-1.2491396) q[1];
sx q[1];
rz(-0.12078665) q[1];
rz(-pi) q[2];
rz(-2.6441908) q[3];
sx q[3];
rz(-0.38555749) q[3];
sx q[3];
rz(1.7413643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5493245) q[2];
sx q[2];
rz(-0.78846875) q[2];
sx q[2];
rz(0.92174021) q[2];
rz(-2.8978469) q[3];
sx q[3];
rz(-1.4344401) q[3];
sx q[3];
rz(2.8488979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.155596) q[0];
sx q[0];
rz(-1.3256185) q[0];
sx q[0];
rz(-2.6407114) q[0];
rz(2.0241375) q[1];
sx q[1];
rz(-1.8084348) q[1];
sx q[1];
rz(-0.31455988) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.019032) q[0];
sx q[0];
rz(-1.9261596) q[0];
sx q[0];
rz(1.2073231) q[0];
rz(-1.6734413) q[2];
sx q[2];
rz(-0.043641239) q[2];
sx q[2];
rz(2.1473644) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.52175888) q[1];
sx q[1];
rz(-1.9653826) q[1];
sx q[1];
rz(0.97418262) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8338449) q[3];
sx q[3];
rz(-0.42908731) q[3];
sx q[3];
rz(-1.3671631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8353117) q[2];
sx q[2];
rz(-1.1113144) q[2];
sx q[2];
rz(1.296898) q[2];
rz(2.6563307) q[3];
sx q[3];
rz(-1.5596215) q[3];
sx q[3];
rz(-2.0280793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6473963) q[0];
sx q[0];
rz(-1.8504471) q[0];
sx q[0];
rz(3.0418292) q[0];
rz(-1.2707204) q[1];
sx q[1];
rz(-2.2427509) q[1];
sx q[1];
rz(-1.3894003) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6806599) q[0];
sx q[0];
rz(-1.4062728) q[0];
sx q[0];
rz(-2.8712832) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9676081) q[2];
sx q[2];
rz(-1.6743039) q[2];
sx q[2];
rz(-3.0826718) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15611616) q[1];
sx q[1];
rz(-2.1735272) q[1];
sx q[1];
rz(-2.3975357) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12677315) q[3];
sx q[3];
rz(-2.354458) q[3];
sx q[3];
rz(1.8925557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.49326593) q[2];
sx q[2];
rz(-0.47097012) q[2];
sx q[2];
rz(-0.89034447) q[2];
rz(-1.9503615) q[3];
sx q[3];
rz(-1.8343364) q[3];
sx q[3];
rz(2.2395535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4297727) q[0];
sx q[0];
rz(-3.0968102) q[0];
sx q[0];
rz(3.1088767) q[0];
rz(-0.50321594) q[1];
sx q[1];
rz(-1.0992173) q[1];
sx q[1];
rz(-0.12282664) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5779992) q[0];
sx q[0];
rz(-1.2482572) q[0];
sx q[0];
rz(2.3911227) q[0];
rz(-pi) q[1];
rz(-1.2314447) q[2];
sx q[2];
rz(-1.095158) q[2];
sx q[2];
rz(-0.95107691) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87007755) q[1];
sx q[1];
rz(-2.6598499) q[1];
sx q[1];
rz(-1.7863356) q[1];
rz(2.8005573) q[3];
sx q[3];
rz(-0.37217316) q[3];
sx q[3];
rz(0.91590524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.34933019) q[2];
sx q[2];
rz(-0.96668875) q[2];
sx q[2];
rz(1.1518504) q[2];
rz(-0.88519111) q[3];
sx q[3];
rz(-1.3693634) q[3];
sx q[3];
rz(-1.7159897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18871466) q[0];
sx q[0];
rz(-2.2080244) q[0];
sx q[0];
rz(-0.92920148) q[0];
rz(0.62905351) q[1];
sx q[1];
rz(-1.7864497) q[1];
sx q[1];
rz(0.74310511) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4152195) q[0];
sx q[0];
rz(-0.45436828) q[0];
sx q[0];
rz(1.3903862) q[0];
rz(2.9628108) q[2];
sx q[2];
rz(-2.6465073) q[2];
sx q[2];
rz(-0.21735969) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3071909) q[1];
sx q[1];
rz(-1.9842923) q[1];
sx q[1];
rz(-1.6430699) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1984897) q[3];
sx q[3];
rz(-2.4022958) q[3];
sx q[3];
rz(-0.044667808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5762081) q[2];
sx q[2];
rz(-2.3614063) q[2];
sx q[2];
rz(-2.4156477) q[2];
rz(-0.37096008) q[3];
sx q[3];
rz(-1.5981263) q[3];
sx q[3];
rz(-0.36271873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87042701) q[0];
sx q[0];
rz(-2.3112264) q[0];
sx q[0];
rz(-2.5776432) q[0];
rz(-1.9922527) q[1];
sx q[1];
rz(-2.1795858) q[1];
sx q[1];
rz(-0.74251485) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19631736) q[0];
sx q[0];
rz(-0.68500297) q[0];
sx q[0];
rz(-2.3514868) q[0];
rz(-pi) q[1];
rz(2.5385802) q[2];
sx q[2];
rz(-2.6508287) q[2];
sx q[2];
rz(0.025207504) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5919784) q[1];
sx q[1];
rz(-1.6149628) q[1];
sx q[1];
rz(0.43818922) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56470676) q[3];
sx q[3];
rz(-0.93834025) q[3];
sx q[3];
rz(-0.45470995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5138862) q[2];
sx q[2];
rz(-1.0739001) q[2];
sx q[2];
rz(2.7970496) q[2];
rz(-2.6978317) q[3];
sx q[3];
rz(-1.0756451) q[3];
sx q[3];
rz(1.3226604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3798856) q[0];
sx q[0];
rz(-0.2825309) q[0];
sx q[0];
rz(2.479582) q[0];
rz(0.33540353) q[1];
sx q[1];
rz(-1.4619724) q[1];
sx q[1];
rz(2.2047156) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7279908) q[0];
sx q[0];
rz(-1.9303891) q[0];
sx q[0];
rz(2.2171564) q[0];
rz(3.1179948) q[2];
sx q[2];
rz(-1.3006217) q[2];
sx q[2];
rz(0.23575704) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.42439991) q[1];
sx q[1];
rz(-2.6584627) q[1];
sx q[1];
rz(-0.65349726) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2912737) q[3];
sx q[3];
rz(-0.26278824) q[3];
sx q[3];
rz(2.9966054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.718049) q[2];
sx q[2];
rz(-2.9073538) q[2];
sx q[2];
rz(-0.49918175) q[2];
rz(1.4623803) q[3];
sx q[3];
rz(-1.9951818) q[3];
sx q[3];
rz(-1.4596938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.725631) q[0];
sx q[0];
rz(-1.220663) q[0];
sx q[0];
rz(2.2882373) q[0];
rz(0.60494963) q[1];
sx q[1];
rz(-1.9974983) q[1];
sx q[1];
rz(-2.6430184) q[1];
rz(-1.7164604) q[2];
sx q[2];
rz(-1.002641) q[2];
sx q[2];
rz(-2.6517131) q[2];
rz(-1.3868757) q[3];
sx q[3];
rz(-1.0763604) q[3];
sx q[3];
rz(3.0767783) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
