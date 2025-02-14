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
rz(5.8402099) q[0];
sx q[0];
rz(9.0415434) q[0];
rz(-1.2869599) q[1];
sx q[1];
rz(1.9998963) q[1];
sx q[1];
rz(13.286446) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38720643) q[0];
sx q[0];
rz(-1.1875679) q[0];
sx q[0];
rz(-2.711074) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.371648) q[2];
sx q[2];
rz(-1.3825644) q[2];
sx q[2];
rz(-3.0668099) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.076700828) q[1];
sx q[1];
rz(-1.2982839) q[1];
sx q[1];
rz(-2.7530159) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96020697) q[3];
sx q[3];
rz(-1.7060602) q[3];
sx q[3];
rz(1.8253411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8522475) q[2];
sx q[2];
rz(-0.87163681) q[2];
sx q[2];
rz(2.3343425) q[2];
rz(1.4489669) q[3];
sx q[3];
rz(-2.2117386) q[3];
sx q[3];
rz(-2.4248031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0221231) q[0];
sx q[0];
rz(-1.2998281) q[0];
sx q[0];
rz(1.8362554) q[0];
rz(-0.7430101) q[1];
sx q[1];
rz(-2.4842255) q[1];
sx q[1];
rz(1.2974799) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9564222) q[0];
sx q[0];
rz(-1.3885986) q[0];
sx q[0];
rz(-1.9745898) q[0];
rz(-pi) q[1];
rz(-2.6941239) q[2];
sx q[2];
rz(-2.0966242) q[2];
sx q[2];
rz(1.4065789) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4549437) q[1];
sx q[1];
rz(-1.8669858) q[1];
sx q[1];
rz(-2.185883) q[1];
rz(-pi) q[2];
rz(1.4543919) q[3];
sx q[3];
rz(-1.9299514) q[3];
sx q[3];
rz(-2.9352279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.3798736) q[2];
sx q[2];
rz(-0.97731176) q[2];
sx q[2];
rz(-1.0404111) q[2];
rz(-0.31446332) q[3];
sx q[3];
rz(-1.9324666) q[3];
sx q[3];
rz(1.8837121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33309015) q[0];
sx q[0];
rz(-0.98121488) q[0];
sx q[0];
rz(1.3370978) q[0];
rz(0.62726504) q[1];
sx q[1];
rz(-1.7341055) q[1];
sx q[1];
rz(-1.6076535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21625265) q[0];
sx q[0];
rz(-2.0006764) q[0];
sx q[0];
rz(2.2735808) q[0];
rz(2.5467403) q[2];
sx q[2];
rz(-0.50083465) q[2];
sx q[2];
rz(-2.5822617) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9107289) q[1];
sx q[1];
rz(-1.7553991) q[1];
sx q[1];
rz(2.6631825) q[1];
x q[2];
rz(-1.1008784) q[3];
sx q[3];
rz(-0.4199314) q[3];
sx q[3];
rz(2.7222999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4072121) q[2];
sx q[2];
rz(-2.247159) q[2];
sx q[2];
rz(1.2073995) q[2];
rz(2.1842128) q[3];
sx q[3];
rz(-1.2355685) q[3];
sx q[3];
rz(-0.68937075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44416881) q[0];
sx q[0];
rz(-0.95132315) q[0];
sx q[0];
rz(3.0498258) q[0];
rz(-2.9755039) q[1];
sx q[1];
rz(-1.2867915) q[1];
sx q[1];
rz(-1.0413569) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4252151) q[0];
sx q[0];
rz(-1.5272523) q[0];
sx q[0];
rz(-2.3460044) q[0];
x q[1];
rz(-0.78775405) q[2];
sx q[2];
rz(-1.279289) q[2];
sx q[2];
rz(2.0306092) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2745121) q[1];
sx q[1];
rz(-1.5364944) q[1];
sx q[1];
rz(2.0438036) q[1];
rz(2.5934173) q[3];
sx q[3];
rz(-0.16528567) q[3];
sx q[3];
rz(-2.838244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.912821) q[2];
sx q[2];
rz(-2.3232338) q[2];
sx q[2];
rz(1.6010326) q[2];
rz(1.0176963) q[3];
sx q[3];
rz(-1.8120268) q[3];
sx q[3];
rz(1.8890007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8215264) q[0];
sx q[0];
rz(-1.7463266) q[0];
sx q[0];
rz(0.23859247) q[0];
rz(0.78856167) q[1];
sx q[1];
rz(-2.8193654) q[1];
sx q[1];
rz(0.90071789) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.032388587) q[0];
sx q[0];
rz(-0.51133832) q[0];
sx q[0];
rz(-0.5503668) q[0];
x q[1];
rz(-0.094036799) q[2];
sx q[2];
rz(-1.9584736) q[2];
sx q[2];
rz(-2.2212818) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2701931) q[1];
sx q[1];
rz(-2.7084623) q[1];
sx q[1];
rz(0.30850839) q[1];
x q[2];
rz(2.2969518) q[3];
sx q[3];
rz(-2.1535404) q[3];
sx q[3];
rz(-2.304519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0476524) q[2];
sx q[2];
rz(-1.36146) q[2];
sx q[2];
rz(1.368604) q[2];
rz(-2.7641344) q[3];
sx q[3];
rz(-0.82337514) q[3];
sx q[3];
rz(1.1641758) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0193598) q[0];
sx q[0];
rz(-2.7074809) q[0];
sx q[0];
rz(-2.9121616) q[0];
rz(-0.58652985) q[1];
sx q[1];
rz(-2.6507381) q[1];
sx q[1];
rz(1.6612926) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0600046) q[0];
sx q[0];
rz(-1.9236653) q[0];
sx q[0];
rz(2.2553483) q[0];
rz(1.3103799) q[2];
sx q[2];
rz(-2.1245109) q[2];
sx q[2];
rz(1.2147724) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.28162128) q[1];
sx q[1];
rz(-0.23275092) q[1];
sx q[1];
rz(-2.908246) q[1];
x q[2];
rz(1.839572) q[3];
sx q[3];
rz(-2.44938) q[3];
sx q[3];
rz(2.6609294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71969879) q[2];
sx q[2];
rz(-2.5514558) q[2];
sx q[2];
rz(1.7936919) q[2];
rz(-0.19139309) q[3];
sx q[3];
rz(-0.32262155) q[3];
sx q[3];
rz(-1.2389244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.557068) q[0];
sx q[0];
rz(-1.5140336) q[0];
sx q[0];
rz(2.6051482) q[0];
rz(-0.084511936) q[1];
sx q[1];
rz(-2.5630496) q[1];
sx q[1];
rz(1.9523778) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7285648) q[0];
sx q[0];
rz(-2.4997093) q[0];
sx q[0];
rz(1.311614) q[0];
rz(-1.9574478) q[2];
sx q[2];
rz(-2.2154567) q[2];
sx q[2];
rz(1.8895011) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.099965) q[1];
sx q[1];
rz(-2.065127) q[1];
sx q[1];
rz(0.35677856) q[1];
rz(-2.0310845) q[3];
sx q[3];
rz(-1.7841859) q[3];
sx q[3];
rz(-1.302857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1875923) q[2];
sx q[2];
rz(-1.8686998) q[2];
sx q[2];
rz(-0.083258955) q[2];
rz(-1.2620874) q[3];
sx q[3];
rz(-0.83674651) q[3];
sx q[3];
rz(2.4038103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5543723) q[0];
sx q[0];
rz(-1.502259) q[0];
sx q[0];
rz(-0.7837522) q[0];
rz(1.4844249) q[1];
sx q[1];
rz(-0.52394167) q[1];
sx q[1];
rz(0.2934244) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3104738) q[0];
sx q[0];
rz(-1.9378512) q[0];
sx q[0];
rz(-1.0793174) q[0];
rz(0.22343724) q[2];
sx q[2];
rz(-0.91470892) q[2];
sx q[2];
rz(0.81435299) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30006584) q[1];
sx q[1];
rz(-1.4264538) q[1];
sx q[1];
rz(1.1529902) q[1];
x q[2];
rz(-2.226172) q[3];
sx q[3];
rz(-1.7599981) q[3];
sx q[3];
rz(-2.9331741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7744814) q[2];
sx q[2];
rz(-1.4159055) q[2];
sx q[2];
rz(1.6488546) q[2];
rz(-1.9914918) q[3];
sx q[3];
rz(-2.8425358) q[3];
sx q[3];
rz(-2.2404631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25006008) q[0];
sx q[0];
rz(-1.0279259) q[0];
sx q[0];
rz(0.0082536396) q[0];
rz(-0.0025657733) q[1];
sx q[1];
rz(-2.6305514) q[1];
sx q[1];
rz(-1.0645688) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7637007) q[0];
sx q[0];
rz(-2.3449909) q[0];
sx q[0];
rz(-2.5515351) q[0];
rz(0.0081811706) q[2];
sx q[2];
rz(-1.3213833) q[2];
sx q[2];
rz(-1.3117664) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9843205) q[1];
sx q[1];
rz(-1.1295348) q[1];
sx q[1];
rz(1.0024789) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5523334) q[3];
sx q[3];
rz(-1.4472826) q[3];
sx q[3];
rz(0.46621399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0801487) q[2];
sx q[2];
rz(-1.4465569) q[2];
sx q[2];
rz(-2.9385369) q[2];
rz(-1.2114581) q[3];
sx q[3];
rz(-1.0695499) q[3];
sx q[3];
rz(-0.13574654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1900629) q[0];
sx q[0];
rz(-1.1030581) q[0];
sx q[0];
rz(0.67361012) q[0];
rz(0.31189648) q[1];
sx q[1];
rz(-1.936828) q[1];
sx q[1];
rz(-1.9685251) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6658026) q[0];
sx q[0];
rz(-0.6190933) q[0];
sx q[0];
rz(1.0503011) q[0];
rz(-pi) q[1];
rz(-0.756097) q[2];
sx q[2];
rz(-1.399516) q[2];
sx q[2];
rz(2.0012228) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7562779) q[1];
sx q[1];
rz(-2.7201001) q[1];
sx q[1];
rz(1.384686) q[1];
x q[2];
rz(-1.0836175) q[3];
sx q[3];
rz(-0.77278944) q[3];
sx q[3];
rz(-0.98675283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6761026) q[2];
sx q[2];
rz(-2.5238621) q[2];
sx q[2];
rz(-0.22107302) q[2];
rz(-2.5681833) q[3];
sx q[3];
rz(-0.87380496) q[3];
sx q[3];
rz(-2.0360086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7318646) q[0];
sx q[0];
rz(-1.8479713) q[0];
sx q[0];
rz(0.027298409) q[0];
rz(-1.191054) q[1];
sx q[1];
rz(-1.4203032) q[1];
sx q[1];
rz(1.3994093) q[1];
rz(0.22573698) q[2];
sx q[2];
rz(-1.0229243) q[2];
sx q[2];
rz(-1.682874) q[2];
rz(-2.1089907) q[3];
sx q[3];
rz(-2.474244) q[3];
sx q[3];
rz(0.49082951) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
