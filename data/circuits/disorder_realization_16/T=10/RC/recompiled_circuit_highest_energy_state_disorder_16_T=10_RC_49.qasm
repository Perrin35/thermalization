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
rz(-1.1666522) q[0];
sx q[0];
rz(-0.44297543) q[0];
sx q[0];
rz(-0.38323453) q[0];
rz(-1.2869599) q[1];
sx q[1];
rz(1.9998963) q[1];
sx q[1];
rz(13.286446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1280649) q[0];
sx q[0];
rz(-1.9682471) q[0];
sx q[0];
rz(-1.1532408) q[0];
rz(-pi) q[1];
rz(2.949657) q[2];
sx q[2];
rz(-1.3752116) q[2];
sx q[2];
rz(1.6078313) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6038415) q[1];
sx q[1];
rz(-1.1972812) q[1];
sx q[1];
rz(-1.2775263) q[1];
rz(-pi) q[2];
rz(-0.96020697) q[3];
sx q[3];
rz(-1.4355324) q[3];
sx q[3];
rz(1.3162515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8522475) q[2];
sx q[2];
rz(-2.2699558) q[2];
sx q[2];
rz(-2.3343425) q[2];
rz(-1.4489669) q[3];
sx q[3];
rz(-0.92985409) q[3];
sx q[3];
rz(-2.4248031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0221231) q[0];
sx q[0];
rz(-1.8417646) q[0];
sx q[0];
rz(1.8362554) q[0];
rz(-0.7430101) q[1];
sx q[1];
rz(-0.65736714) q[1];
sx q[1];
rz(-1.2974799) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8332307) q[0];
sx q[0];
rz(-1.9675281) q[0];
sx q[0];
rz(-2.9438567) q[0];
rz(-pi) q[1];
rz(-2.6941239) q[2];
sx q[2];
rz(-1.0449685) q[2];
sx q[2];
rz(1.7350137) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.686649) q[1];
sx q[1];
rz(-1.8669858) q[1];
sx q[1];
rz(0.95570968) q[1];
x q[2];
rz(-2.8415806) q[3];
sx q[3];
rz(-2.7648246) q[3];
sx q[3];
rz(3.026791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7617191) q[2];
sx q[2];
rz(-0.97731176) q[2];
sx q[2];
rz(2.1011815) q[2];
rz(-0.31446332) q[3];
sx q[3];
rz(-1.9324666) q[3];
sx q[3];
rz(-1.2578806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8085025) q[0];
sx q[0];
rz(-0.98121488) q[0];
sx q[0];
rz(1.3370978) q[0];
rz(-0.62726504) q[1];
sx q[1];
rz(-1.4074872) q[1];
sx q[1];
rz(1.5339392) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8117789) q[0];
sx q[0];
rz(-2.3372991) q[0];
sx q[0];
rz(-0.95383184) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4256469) q[2];
sx q[2];
rz(-1.2983649) q[2];
sx q[2];
rz(-2.6657588) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1418504) q[1];
sx q[1];
rz(-2.6313884) q[1];
sx q[1];
rz(-0.38533089) q[1];
rz(-1.1919349) q[3];
sx q[3];
rz(-1.7564723) q[3];
sx q[3];
rz(2.4243086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73438054) q[2];
sx q[2];
rz(-0.89443365) q[2];
sx q[2];
rz(1.2073995) q[2];
rz(-2.1842128) q[3];
sx q[3];
rz(-1.9060241) q[3];
sx q[3];
rz(-0.68937075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44416881) q[0];
sx q[0];
rz(-0.95132315) q[0];
sx q[0];
rz(0.091766894) q[0];
rz(0.16608876) q[1];
sx q[1];
rz(-1.8548012) q[1];
sx q[1];
rz(1.0413569) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7163775) q[0];
sx q[0];
rz(-1.6143403) q[0];
sx q[0];
rz(-0.79558827) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40046863) q[2];
sx q[2];
rz(-2.31268) q[2];
sx q[2];
rz(0.18120719) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68616679) q[1];
sx q[1];
rz(-1.09809) q[1];
sx q[1];
rz(0.03852877) q[1];
rz(-1.4840858) q[3];
sx q[3];
rz(-1.711688) q[3];
sx q[3];
rz(2.2839462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.2287717) q[2];
sx q[2];
rz(-2.3232338) q[2];
sx q[2];
rz(-1.6010326) q[2];
rz(2.1238964) q[3];
sx q[3];
rz(-1.3295659) q[3];
sx q[3];
rz(1.8890007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3200662) q[0];
sx q[0];
rz(-1.3952661) q[0];
sx q[0];
rz(0.23859247) q[0];
rz(-2.353031) q[1];
sx q[1];
rz(-0.32222727) q[1];
sx q[1];
rz(-0.90071789) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.032388587) q[0];
sx q[0];
rz(-0.51133832) q[0];
sx q[0];
rz(-2.5912259) q[0];
rz(-1.7968171) q[2];
sx q[2];
rz(-0.39835762) q[2];
sx q[2];
rz(-2.4657754) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93253329) q[1];
sx q[1];
rz(-1.9822023) q[1];
sx q[1];
rz(1.4313015) q[1];
rz(-0.84464083) q[3];
sx q[3];
rz(-0.9880522) q[3];
sx q[3];
rz(2.304519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0476524) q[2];
sx q[2];
rz(-1.7801327) q[2];
sx q[2];
rz(1.7729887) q[2];
rz(2.7641344) q[3];
sx q[3];
rz(-0.82337514) q[3];
sx q[3];
rz(-1.1641758) q[3];
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
rz(-2.1222328) q[0];
sx q[0];
rz(-2.7074809) q[0];
sx q[0];
rz(-2.9121616) q[0];
rz(2.5550628) q[1];
sx q[1];
rz(-0.49085453) q[1];
sx q[1];
rz(-1.6612926) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2302814) q[0];
sx q[0];
rz(-2.3846941) q[0];
sx q[0];
rz(-2.0981936) q[0];
rz(1.3103799) q[2];
sx q[2];
rz(-2.1245109) q[2];
sx q[2];
rz(-1.9268203) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0796894) q[1];
sx q[1];
rz(-1.6241571) q[1];
sx q[1];
rz(-2.9149331) q[1];
rz(-pi) q[2];
x q[2];
rz(1.839572) q[3];
sx q[3];
rz(-2.44938) q[3];
sx q[3];
rz(-0.48066329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4218939) q[2];
sx q[2];
rz(-0.59013683) q[2];
sx q[2];
rz(-1.7936919) q[2];
rz(2.9501996) q[3];
sx q[3];
rz(-0.32262155) q[3];
sx q[3];
rz(-1.2389244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58452463) q[0];
sx q[0];
rz(-1.5140336) q[0];
sx q[0];
rz(2.6051482) q[0];
rz(3.0570807) q[1];
sx q[1];
rz(-0.5785431) q[1];
sx q[1];
rz(-1.9523778) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41302785) q[0];
sx q[0];
rz(-2.4997093) q[0];
sx q[0];
rz(1.8299787) q[0];
rz(0.68185735) q[2];
sx q[2];
rz(-1.2646265) q[2];
sx q[2];
rz(-0.078753565) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0416276) q[1];
sx q[1];
rz(-1.0764657) q[1];
sx q[1];
rz(2.7848141) q[1];
rz(-pi) q[2];
rz(-2.9042894) q[3];
sx q[3];
rz(-1.1217228) q[3];
sx q[3];
rz(-0.16332353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1875923) q[2];
sx q[2];
rz(-1.2728929) q[2];
sx q[2];
rz(0.083258955) q[2];
rz(1.2620874) q[3];
sx q[3];
rz(-2.3048461) q[3];
sx q[3];
rz(2.4038103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5543723) q[0];
sx q[0];
rz(-1.6393336) q[0];
sx q[0];
rz(0.7837522) q[0];
rz(1.6571677) q[1];
sx q[1];
rz(-2.617651) q[1];
sx q[1];
rz(0.2934244) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2121289) q[0];
sx q[0];
rz(-1.1146472) q[0];
sx q[0];
rz(-0.41123234) q[0];
x q[1];
rz(0.90250315) q[2];
sx q[2];
rz(-1.74729) q[2];
sx q[2];
rz(0.89418156) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95739505) q[1];
sx q[1];
rz(-0.44064545) q[1];
sx q[1];
rz(-1.9147826) q[1];
rz(-pi) q[2];
rz(2.226172) q[3];
sx q[3];
rz(-1.3815945) q[3];
sx q[3];
rz(0.20841852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.36711127) q[2];
sx q[2];
rz(-1.7256871) q[2];
sx q[2];
rz(1.492738) q[2];
rz(-1.9914918) q[3];
sx q[3];
rz(-2.8425358) q[3];
sx q[3];
rz(-2.2404631) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8915326) q[0];
sx q[0];
rz(-2.1136668) q[0];
sx q[0];
rz(-3.133339) q[0];
rz(3.1390269) q[1];
sx q[1];
rz(-2.6305514) q[1];
sx q[1];
rz(-1.0645688) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9998523) q[0];
sx q[0];
rz(-0.93466648) q[0];
sx q[0];
rz(-1.0534691) q[0];
rz(-pi) q[1];
rz(-0.0081811706) q[2];
sx q[2];
rz(-1.8202094) q[2];
sx q[2];
rz(1.8298263) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1572722) q[1];
sx q[1];
rz(-1.1295348) q[1];
sx q[1];
rz(-2.1391137) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4225619) q[3];
sx q[3];
rz(-2.1549716) q[3];
sx q[3];
rz(2.119182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0801487) q[2];
sx q[2];
rz(-1.6950357) q[2];
sx q[2];
rz(-0.20305571) q[2];
rz(-1.2114581) q[3];
sx q[3];
rz(-1.0695499) q[3];
sx q[3];
rz(3.0058461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-1.9515297) q[0];
sx q[0];
rz(-1.1030581) q[0];
sx q[0];
rz(2.4679825) q[0];
rz(-2.8296962) q[1];
sx q[1];
rz(-1.936828) q[1];
sx q[1];
rz(-1.9685251) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53176994) q[0];
sx q[0];
rz(-1.2780452) q[0];
sx q[0];
rz(-1.0171153) q[0];
rz(-pi) q[1];
rz(-2.8946213) q[2];
sx q[2];
rz(-2.3701043) q[2];
sx q[2];
rz(-2.889973) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0153326) q[1];
sx q[1];
rz(-1.646572) q[1];
sx q[1];
rz(-1.9858236) q[1];
x q[2];
rz(-0.85957771) q[3];
sx q[3];
rz(-1.2378581) q[3];
sx q[3];
rz(-2.9200819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6761026) q[2];
sx q[2];
rz(-2.5238621) q[2];
sx q[2];
rz(-0.22107302) q[2];
rz(0.57340932) q[3];
sx q[3];
rz(-2.2677877) q[3];
sx q[3];
rz(2.0360086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7318646) q[0];
sx q[0];
rz(-1.8479713) q[0];
sx q[0];
rz(0.027298409) q[0];
rz(1.191054) q[1];
sx q[1];
rz(-1.7212894) q[1];
sx q[1];
rz(-1.7421834) q[1];
rz(1.9223735) q[2];
sx q[2];
rz(-2.5534591) q[2];
sx q[2];
rz(-2.0981325) q[2];
rz(-0.97597573) q[3];
sx q[3];
rz(-1.2479758) q[3];
sx q[3];
rz(1.6231619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
