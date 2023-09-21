OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1749984) q[0];
sx q[0];
rz(2.7881665) q[0];
sx q[0];
rz(8.3600144) q[0];
rz(0.79617533) q[1];
sx q[1];
rz(-1.9328971) q[1];
sx q[1];
rz(-2.605521) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85815) q[0];
sx q[0];
rz(-1.3387696) q[0];
sx q[0];
rz(2.0767077) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0739003) q[2];
sx q[2];
rz(-2.2798685) q[2];
sx q[2];
rz(-0.46082218) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88692611) q[1];
sx q[1];
rz(-0.98985043) q[1];
sx q[1];
rz(1.7960153) q[1];
x q[2];
rz(0.73379559) q[3];
sx q[3];
rz(-1.9473837) q[3];
sx q[3];
rz(0.22613444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7321695) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(-2.1146963) q[2];
rz(-2.8895767) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(-0.86565971) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058763) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(-1.3902364) q[0];
rz(0.83579666) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(-0.70835152) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81472155) q[0];
sx q[0];
rz(-0.46778361) q[0];
sx q[0];
rz(-0.65713709) q[0];
rz(-2.9950036) q[2];
sx q[2];
rz(-0.97448889) q[2];
sx q[2];
rz(1.6080315) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2409776) q[1];
sx q[1];
rz(-2.2574189) q[1];
sx q[1];
rz(2.2254506) q[1];
x q[2];
rz(2.0833011) q[3];
sx q[3];
rz(-1.1304434) q[3];
sx q[3];
rz(-3.1194221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11671242) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(2.8498245) q[2];
rz(3.0388888) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(-1.5244938) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8125238) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(2.8318751) q[0];
rz(-1.4801056) q[1];
sx q[1];
rz(-1.8146351) q[1];
sx q[1];
rz(2.5699239) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0535307) q[0];
sx q[0];
rz(-0.54170875) q[0];
sx q[0];
rz(2.9689229) q[0];
rz(-pi) q[1];
rz(1.13582) q[2];
sx q[2];
rz(-1.2425213) q[2];
sx q[2];
rz(0.091718397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7844312) q[1];
sx q[1];
rz(-2.2964587) q[1];
sx q[1];
rz(1.3933338) q[1];
x q[2];
rz(0.27922697) q[3];
sx q[3];
rz(-0.46776566) q[3];
sx q[3];
rz(1.8750909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16033515) q[2];
sx q[2];
rz(-2.2312639) q[2];
sx q[2];
rz(-2.4327915) q[2];
rz(2.839084) q[3];
sx q[3];
rz(-1.442028) q[3];
sx q[3];
rz(-1.1423053) q[3];
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
rz(-0.70808327) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(2.3420912) q[0];
rz(3.0917621) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(-2.9611011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7364396) q[0];
sx q[0];
rz(-2.7442051) q[0];
sx q[0];
rz(2.9050164) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7454342) q[2];
sx q[2];
rz(-1.9348113) q[2];
sx q[2];
rz(-2.3166235) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0536249) q[1];
sx q[1];
rz(-1.2521724) q[1];
sx q[1];
rz(1.6941403) q[1];
rz(-1.5482076) q[3];
sx q[3];
rz(-0.49593192) q[3];
sx q[3];
rz(-3.0043234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0522456) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(-1.583741) q[2];
rz(-1.3752939) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-0.93262514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76064008) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(0.98651648) q[0];
rz(1.1622693) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(-0.25156897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7154327) q[0];
sx q[0];
rz(-2.4315205) q[0];
sx q[0];
rz(-1.8415113) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6576084) q[2];
sx q[2];
rz(-2.0292536) q[2];
sx q[2];
rz(0.77755962) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0136452) q[1];
sx q[1];
rz(-0.19731678) q[1];
sx q[1];
rz(3.0017588) q[1];
x q[2];
rz(-2.6336446) q[3];
sx q[3];
rz(-2.2077999) q[3];
sx q[3];
rz(0.38364832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(-1.0305369) q[2];
rz(1.7306227) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(-0.63123909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3866766) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(2.2914698) q[0];
rz(1.9592346) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(-2.7485671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5315006) q[0];
sx q[0];
rz(-1.8787662) q[0];
sx q[0];
rz(-2.969335) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23586258) q[2];
sx q[2];
rz(-0.84025192) q[2];
sx q[2];
rz(-2.0116531) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1054354) q[1];
sx q[1];
rz(-1.565355) q[1];
sx q[1];
rz(2.2219707) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29720184) q[3];
sx q[3];
rz(-1.3085877) q[3];
sx q[3];
rz(-2.8924243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(0.97314107) q[2];
rz(2.1940103) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(-1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.4457552) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(-1.6181035) q[0];
rz(2.7667926) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(0.0059676776) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25516605) q[0];
sx q[0];
rz(-0.92029858) q[0];
sx q[0];
rz(-2.7889473) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75255021) q[2];
sx q[2];
rz(-1.2298905) q[2];
sx q[2];
rz(-0.37170751) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.70049268) q[1];
sx q[1];
rz(-1.0171434) q[1];
sx q[1];
rz(-3.0667552) q[1];
x q[2];
rz(0.11984101) q[3];
sx q[3];
rz(-0.37444515) q[3];
sx q[3];
rz(-0.10711174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.57198793) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(2.9515284) q[2];
rz(2.7251785) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(0.73474187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1104601) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(-0.53623143) q[0];
rz(2.2082632) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(-0.94820625) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1873916) q[0];
sx q[0];
rz(-1.3905977) q[0];
sx q[0];
rz(-0.76822922) q[0];
rz(-1.0529222) q[2];
sx q[2];
rz(-1.1793943) q[2];
sx q[2];
rz(-0.97757593) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.58029786) q[1];
sx q[1];
rz(-1.6333688) q[1];
sx q[1];
rz(-0.22217521) q[1];
x q[2];
rz(-0.35685278) q[3];
sx q[3];
rz(-1.4992504) q[3];
sx q[3];
rz(-2.1944291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.640921) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(-0.24027696) q[2];
rz(2.7622973) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(-1.3482288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.76294476) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(-2.8906524) q[0];
rz(-2.8885686) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(-0.35266638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6610722) q[0];
sx q[0];
rz(-1.4561704) q[0];
sx q[0];
rz(-1.4474439) q[0];
x q[1];
rz(-0.30013957) q[2];
sx q[2];
rz(-1.7087414) q[2];
sx q[2];
rz(2.2696242) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.36665146) q[1];
sx q[1];
rz(-2.1244441) q[1];
sx q[1];
rz(-1.8371546) q[1];
rz(-pi) q[2];
rz(-0.50036939) q[3];
sx q[3];
rz(-2.2420068) q[3];
sx q[3];
rz(-2.3994115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76715604) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(-2.4582668) q[2];
rz(1.654489) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(-1.1504415) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3764573) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(-1.8413683) q[0];
rz(-0.70872778) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(-0.93651071) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4721601) q[0];
sx q[0];
rz(-0.59211187) q[0];
sx q[0];
rz(-2.6955312) q[0];
rz(-0.18239393) q[2];
sx q[2];
rz(-1.2904022) q[2];
sx q[2];
rz(-0.80255752) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6259389) q[1];
sx q[1];
rz(-1.7682372) q[1];
sx q[1];
rz(-2.291912) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1533222) q[3];
sx q[3];
rz(-1.3068145) q[3];
sx q[3];
rz(1.9895944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5986754) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(-3.1266406) q[2];
rz(2.9144918) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(-1.2623513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99075714) q[0];
sx q[0];
rz(-2.3322454) q[0];
sx q[0];
rz(-1.1582751) q[0];
rz(-1.5785718) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(1.0583393) q[2];
sx q[2];
rz(-0.28596157) q[2];
sx q[2];
rz(2.7059976) q[2];
rz(-0.43185497) q[3];
sx q[3];
rz(-2.5681744) q[3];
sx q[3];
rz(-1.4011866) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];