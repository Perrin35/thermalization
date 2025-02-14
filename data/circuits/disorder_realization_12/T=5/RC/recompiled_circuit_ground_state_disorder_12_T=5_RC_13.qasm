OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0924858) q[0];
sx q[0];
rz(-1.2993113) q[0];
sx q[0];
rz(3.0957676) q[0];
rz(1.2996281) q[1];
sx q[1];
rz(-1.8027432) q[1];
sx q[1];
rz(1.2531228) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6413413) q[0];
sx q[0];
rz(-1.1317028) q[0];
sx q[0];
rz(2.8706949) q[0];
rz(2.9304977) q[2];
sx q[2];
rz(-2.935754) q[2];
sx q[2];
rz(1.4343651) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7980849) q[1];
sx q[1];
rz(-1.4516593) q[1];
sx q[1];
rz(0.75372523) q[1];
rz(-pi) q[2];
rz(0.0039611749) q[3];
sx q[3];
rz(-1.4433493) q[3];
sx q[3];
rz(-2.733903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8200298) q[2];
sx q[2];
rz(-1.7367312) q[2];
sx q[2];
rz(2.826214) q[2];
rz(0.086221181) q[3];
sx q[3];
rz(-0.092970522) q[3];
sx q[3];
rz(2.2975547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.740199) q[0];
sx q[0];
rz(-1.4305776) q[0];
sx q[0];
rz(-2.8061818) q[0];
rz(-2.7566578) q[1];
sx q[1];
rz(-0.79610577) q[1];
sx q[1];
rz(1.8523432) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9612564) q[0];
sx q[0];
rz(-1.5085876) q[0];
sx q[0];
rz(0.38952413) q[0];
rz(-0.92929594) q[2];
sx q[2];
rz(-1.6469064) q[2];
sx q[2];
rz(-0.92744213) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.030176) q[1];
sx q[1];
rz(-1.9523638) q[1];
sx q[1];
rz(-1.2658582) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13522526) q[3];
sx q[3];
rz(-1.8117419) q[3];
sx q[3];
rz(1.8640445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5424767) q[2];
sx q[2];
rz(-1.7581538) q[2];
sx q[2];
rz(1.6790338) q[2];
rz(-2.2595432) q[3];
sx q[3];
rz(-2.4817011) q[3];
sx q[3];
rz(1.5006458) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3682692) q[0];
sx q[0];
rz(-0.48650807) q[0];
sx q[0];
rz(0.64273709) q[0];
rz(0.87359387) q[1];
sx q[1];
rz(-1.830955) q[1];
sx q[1];
rz(-2.1547623) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9844279) q[0];
sx q[0];
rz(-1.4207977) q[0];
sx q[0];
rz(2.2638129) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0552733) q[2];
sx q[2];
rz(-1.5159801) q[2];
sx q[2];
rz(1.1940741) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1912811) q[1];
sx q[1];
rz(-1.0265959) q[1];
sx q[1];
rz(-0.64375513) q[1];
rz(-pi) q[2];
rz(-0.049792265) q[3];
sx q[3];
rz(-1.9450499) q[3];
sx q[3];
rz(-2.5516537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4387536) q[2];
sx q[2];
rz(-0.025391014) q[2];
sx q[2];
rz(1.0566443) q[2];
rz(1.7993401) q[3];
sx q[3];
rz(-1.6868351) q[3];
sx q[3];
rz(0.87856436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1037647) q[0];
sx q[0];
rz(-2.8570638) q[0];
sx q[0];
rz(-1.1430662) q[0];
rz(-0.68483886) q[1];
sx q[1];
rz(-1.3416483) q[1];
sx q[1];
rz(-2.8917704) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.095465) q[0];
sx q[0];
rz(-2.2860512) q[0];
sx q[0];
rz(-0.94761316) q[0];
rz(-0.5742214) q[2];
sx q[2];
rz(-1.5674233) q[2];
sx q[2];
rz(0.33004883) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4433684) q[1];
sx q[1];
rz(-2.0495546) q[1];
sx q[1];
rz(2.4633292) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8127727) q[3];
sx q[3];
rz(-1.5692838) q[3];
sx q[3];
rz(-2.5082299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7217241) q[2];
sx q[2];
rz(-1.7896174) q[2];
sx q[2];
rz(1.8458337) q[2];
rz(-1.3755679) q[3];
sx q[3];
rz(-1.7121168) q[3];
sx q[3];
rz(3.0425369) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7527723) q[0];
sx q[0];
rz(-1.3714014) q[0];
sx q[0];
rz(-0.8465299) q[0];
rz(-1.1835774) q[1];
sx q[1];
rz(-1.1776244) q[1];
sx q[1];
rz(-1.2394946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1355818) q[0];
sx q[0];
rz(-2.860489) q[0];
sx q[0];
rz(-2.0182039) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34244142) q[2];
sx q[2];
rz(-1.3991809) q[2];
sx q[2];
rz(1.1521074) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2927276) q[1];
sx q[1];
rz(-1.5588817) q[1];
sx q[1];
rz(-1.5148276) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1679513) q[3];
sx q[3];
rz(-1.6203364) q[3];
sx q[3];
rz(1.3508391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5014629) q[2];
sx q[2];
rz(-1.632246) q[2];
sx q[2];
rz(-0.35422361) q[2];
rz(0.7192449) q[3];
sx q[3];
rz(-2.5651599) q[3];
sx q[3];
rz(0.80938068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7648014) q[0];
sx q[0];
rz(-0.23623315) q[0];
sx q[0];
rz(2.7948622) q[0];
rz(2.4193343) q[1];
sx q[1];
rz(-1.5303231) q[1];
sx q[1];
rz(2.9647656) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051561616) q[0];
sx q[0];
rz(-1.430184) q[0];
sx q[0];
rz(1.1052126) q[0];
rz(3.0484285) q[2];
sx q[2];
rz(-2.0679065) q[2];
sx q[2];
rz(2.1456631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74048282) q[1];
sx q[1];
rz(-1.4663457) q[1];
sx q[1];
rz(-2.9354276) q[1];
x q[2];
rz(0.10156722) q[3];
sx q[3];
rz(-1.0103161) q[3];
sx q[3];
rz(1.5965726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56883812) q[2];
sx q[2];
rz(-2.4080031) q[2];
sx q[2];
rz(2.0983762) q[2];
rz(-0.15626945) q[3];
sx q[3];
rz(-1.0633435) q[3];
sx q[3];
rz(0.094303057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909376) q[0];
sx q[0];
rz(-1.0303048) q[0];
sx q[0];
rz(-1.2283196) q[0];
rz(0.43371513) q[1];
sx q[1];
rz(-0.80077306) q[1];
sx q[1];
rz(2.0965651) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50253326) q[0];
sx q[0];
rz(-2.9307588) q[0];
sx q[0];
rz(1.8502214) q[0];
x q[1];
rz(2.3190193) q[2];
sx q[2];
rz(-1.9299091) q[2];
sx q[2];
rz(0.99701478) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.0050538926) q[1];
sx q[1];
rz(-0.95816678) q[1];
sx q[1];
rz(3.0023252) q[1];
rz(-pi) q[2];
rz(-0.52900903) q[3];
sx q[3];
rz(-2.2914791) q[3];
sx q[3];
rz(-1.0320272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9016483) q[2];
sx q[2];
rz(-2.1176691) q[2];
sx q[2];
rz(-2.955692) q[2];
rz(1.2402041) q[3];
sx q[3];
rz(-1.600772) q[3];
sx q[3];
rz(-1.0143636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1321201) q[0];
sx q[0];
rz(-1.2484231) q[0];
sx q[0];
rz(-0.15952071) q[0];
rz(2.3347143) q[1];
sx q[1];
rz(-1.9069549) q[1];
sx q[1];
rz(-0.0085011403) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5497564) q[0];
sx q[0];
rz(-1.4167794) q[0];
sx q[0];
rz(1.9817673) q[0];
rz(-1.6775756) q[2];
sx q[2];
rz(-0.77645436) q[2];
sx q[2];
rz(1.3173033) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2841683) q[1];
sx q[1];
rz(-1.7021952) q[1];
sx q[1];
rz(-2.9144276) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0866585) q[3];
sx q[3];
rz(-1.7268506) q[3];
sx q[3];
rz(-0.54359667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.27935091) q[2];
sx q[2];
rz(-1.3631577) q[2];
sx q[2];
rz(-2.9478493) q[2];
rz(1.5396317) q[3];
sx q[3];
rz(-1.1774457) q[3];
sx q[3];
rz(-2.0120473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9047456) q[0];
sx q[0];
rz(-1.3847677) q[0];
sx q[0];
rz(2.4054476) q[0];
rz(0.88788095) q[1];
sx q[1];
rz(-0.45279756) q[1];
sx q[1];
rz(-0.050051659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56704484) q[0];
sx q[0];
rz(-2.7549681) q[0];
sx q[0];
rz(-0.48166122) q[0];
x q[1];
rz(0.31466575) q[2];
sx q[2];
rz(-1.6648714) q[2];
sx q[2];
rz(-1.5794181) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1117275) q[1];
sx q[1];
rz(-1.3674515) q[1];
sx q[1];
rz(-0.3158658) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1878236) q[3];
sx q[3];
rz(-2.0514384) q[3];
sx q[3];
rz(-2.1565089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2729518) q[2];
sx q[2];
rz(-1.7014039) q[2];
sx q[2];
rz(3.1325373) q[2];
rz(2.791259) q[3];
sx q[3];
rz(-2.83559) q[3];
sx q[3];
rz(0.30456021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.7250799) q[0];
sx q[0];
rz(-1.059499) q[0];
sx q[0];
rz(-2.1685725) q[0];
rz(-1.579772) q[1];
sx q[1];
rz(-0.90619722) q[1];
sx q[1];
rz(0.80950338) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1081297) q[0];
sx q[0];
rz(-1.9282883) q[0];
sx q[0];
rz(0.78127677) q[0];
x q[1];
rz(-1.8000748) q[2];
sx q[2];
rz(-1.7588758) q[2];
sx q[2];
rz(-0.089872472) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0187518) q[1];
sx q[1];
rz(-2.0301308) q[1];
sx q[1];
rz(-0.18045119) q[1];
x q[2];
rz(-2.6736835) q[3];
sx q[3];
rz(-1.0413326) q[3];
sx q[3];
rz(-0.61509815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1805264) q[2];
sx q[2];
rz(-0.82225353) q[2];
sx q[2];
rz(-0.52214885) q[2];
rz(-0.3565878) q[3];
sx q[3];
rz(-0.57727376) q[3];
sx q[3];
rz(-0.063551158) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6575573) q[0];
sx q[0];
rz(-0.35150305) q[0];
sx q[0];
rz(1.2765314) q[0];
rz(-0.36318489) q[1];
sx q[1];
rz(-2.6138432) q[1];
sx q[1];
rz(1.6539727) q[1];
rz(-2.364306) q[2];
sx q[2];
rz(-1.9718134) q[2];
sx q[2];
rz(2.8169463) q[2];
rz(1.894886) q[3];
sx q[3];
rz(-1.0151328) q[3];
sx q[3];
rz(-0.48965164) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
