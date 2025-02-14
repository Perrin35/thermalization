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
rz(-0.18345565) q[0];
sx q[0];
rz(-0.7440716) q[0];
sx q[0];
rz(2.0896572) q[0];
rz(-2.9479041) q[1];
sx q[1];
rz(-2.5084578) q[1];
sx q[1];
rz(2.5685891) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6010701) q[0];
sx q[0];
rz(-2.1946476) q[0];
sx q[0];
rz(-3.0649351) q[0];
x q[1];
rz(-2.6296205) q[2];
sx q[2];
rz(-2.3880771) q[2];
sx q[2];
rz(-1.3775502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6261648) q[1];
sx q[1];
rz(-1.0744796) q[1];
sx q[1];
rz(1.650701) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4450847) q[3];
sx q[3];
rz(-1.4918431) q[3];
sx q[3];
rz(-2.1895308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8594592) q[2];
sx q[2];
rz(-2.8318475) q[2];
sx q[2];
rz(-2.0852883) q[2];
rz(0.022196444) q[3];
sx q[3];
rz(-0.76265097) q[3];
sx q[3];
rz(1.4200042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9098772) q[0];
sx q[0];
rz(-2.660399) q[0];
sx q[0];
rz(-0.43014446) q[0];
rz(-3.0138956) q[1];
sx q[1];
rz(-1.1558497) q[1];
sx q[1];
rz(-1.7040303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6352954) q[0];
sx q[0];
rz(-0.40182913) q[0];
sx q[0];
rz(0.70962972) q[0];
rz(2.2980437) q[2];
sx q[2];
rz(-3.0650716) q[2];
sx q[2];
rz(-0.25329548) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3900657) q[1];
sx q[1];
rz(-1.8987149) q[1];
sx q[1];
rz(-0.71050127) q[1];
rz(-0.35234612) q[3];
sx q[3];
rz(-2.5735613) q[3];
sx q[3];
rz(1.2888792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11640707) q[2];
sx q[2];
rz(-1.4613287) q[2];
sx q[2];
rz(-2.6961668) q[2];
rz(-2.7323501) q[3];
sx q[3];
rz(-1.0990812) q[3];
sx q[3];
rz(2.2621431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5889848) q[0];
sx q[0];
rz(-0.092985066) q[0];
sx q[0];
rz(2.4267922) q[0];
rz(1.0572761) q[1];
sx q[1];
rz(-0.42066586) q[1];
sx q[1];
rz(-0.32726273) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0258758) q[0];
sx q[0];
rz(-2.1109606) q[0];
sx q[0];
rz(-2.9469423) q[0];
x q[1];
rz(-2.1959841) q[2];
sx q[2];
rz(-0.66788061) q[2];
sx q[2];
rz(0.57648522) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5048557) q[1];
sx q[1];
rz(-1.5914006) q[1];
sx q[1];
rz(-1.4192753) q[1];
x q[2];
rz(-2.1932262) q[3];
sx q[3];
rz(-0.41928681) q[3];
sx q[3];
rz(-2.6635827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0032234) q[2];
sx q[2];
rz(-0.97347632) q[2];
sx q[2];
rz(1.5039911) q[2];
rz(1.9231632) q[3];
sx q[3];
rz(-2.7259493) q[3];
sx q[3];
rz(2.159582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0686491) q[0];
sx q[0];
rz(-1.8953841) q[0];
sx q[0];
rz(-2.9856227) q[0];
rz(2.7929557) q[1];
sx q[1];
rz(-0.60395423) q[1];
sx q[1];
rz(1.9085931) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5208682) q[0];
sx q[0];
rz(-1.8308795) q[0];
sx q[0];
rz(3.0164032) q[0];
x q[1];
rz(-1.8470339) q[2];
sx q[2];
rz(-2.4748908) q[2];
sx q[2];
rz(1.897097) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.33267) q[1];
sx q[1];
rz(-2.4415209) q[1];
sx q[1];
rz(-1.1247404) q[1];
rz(-pi) q[2];
rz(3.0844131) q[3];
sx q[3];
rz(-2.8897277) q[3];
sx q[3];
rz(2.4185023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6984581) q[2];
sx q[2];
rz(-3.105574) q[2];
sx q[2];
rz(-1.5114463) q[2];
rz(2.002142) q[3];
sx q[3];
rz(-1.3316493) q[3];
sx q[3];
rz(-0.99036923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.776942) q[0];
sx q[0];
rz(-0.41780892) q[0];
sx q[0];
rz(-1.9885709) q[0];
rz(-2.5979089) q[1];
sx q[1];
rz(-1.8749571) q[1];
sx q[1];
rz(0.26184729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1516929) q[0];
sx q[0];
rz(-1.4873056) q[0];
sx q[0];
rz(0.033647353) q[0];
rz(-pi) q[1];
rz(-2.8072281) q[2];
sx q[2];
rz(-1.0131256) q[2];
sx q[2];
rz(1.1209436) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5763229) q[1];
sx q[1];
rz(-2.5082631) q[1];
sx q[1];
rz(-2.759658) q[1];
rz(-2.9363067) q[3];
sx q[3];
rz(-0.79278273) q[3];
sx q[3];
rz(1.0012116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5220945) q[2];
sx q[2];
rz(-2.5358989) q[2];
sx q[2];
rz(-1.4052793) q[2];
rz(2.3960579) q[3];
sx q[3];
rz(-1.2657974) q[3];
sx q[3];
rz(-2.1982927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.140542) q[0];
sx q[0];
rz(-3.0241835) q[0];
sx q[0];
rz(0.35181272) q[0];
rz(0.44031269) q[1];
sx q[1];
rz(-1.3654717) q[1];
sx q[1];
rz(0.17131677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650314) q[0];
sx q[0];
rz(-0.91463551) q[0];
sx q[0];
rz(-1.879346) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2662884) q[2];
sx q[2];
rz(-0.84104462) q[2];
sx q[2];
rz(2.8685121) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3503987) q[1];
sx q[1];
rz(-1.5499252) q[1];
sx q[1];
rz(-0.69857614) q[1];
x q[2];
rz(0.76254179) q[3];
sx q[3];
rz(-2.2665215) q[3];
sx q[3];
rz(2.1660337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0773086) q[2];
sx q[2];
rz(-1.7694387) q[2];
sx q[2];
rz(-0.94179955) q[2];
rz(1.1549548) q[3];
sx q[3];
rz(-0.20808163) q[3];
sx q[3];
rz(-1.340516) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6963541) q[0];
sx q[0];
rz(-1.1055163) q[0];
sx q[0];
rz(0.0048986991) q[0];
rz(-0.44662961) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(-2.4536224) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1991581) q[0];
sx q[0];
rz(-0.73478886) q[0];
sx q[0];
rz(-0.25811974) q[0];
rz(-pi) q[1];
rz(0.46779386) q[2];
sx q[2];
rz(-2.4681598) q[2];
sx q[2];
rz(-1.837567) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0258642) q[1];
sx q[1];
rz(-2.0021571) q[1];
sx q[1];
rz(-0.23484767) q[1];
x q[2];
rz(2.2673685) q[3];
sx q[3];
rz(-1.0099353) q[3];
sx q[3];
rz(-2.0340597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61116162) q[2];
sx q[2];
rz(-2.7335584) q[2];
sx q[2];
rz(-0.92998695) q[2];
rz(1.9200578) q[3];
sx q[3];
rz(-2.0285716) q[3];
sx q[3];
rz(-0.59337029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68194836) q[0];
sx q[0];
rz(-1.3197897) q[0];
sx q[0];
rz(3.0514858) q[0];
rz(-1.7168761) q[1];
sx q[1];
rz(-2.5017891) q[1];
sx q[1];
rz(-0.43509126) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1662707) q[0];
sx q[0];
rz(-1.4251801) q[0];
sx q[0];
rz(-2.6377489) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66799156) q[2];
sx q[2];
rz(-1.4282303) q[2];
sx q[2];
rz(0.1145471) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4415978) q[1];
sx q[1];
rz(-2.1837103) q[1];
sx q[1];
rz(-2.9843778) q[1];
x q[2];
rz(1.5716721) q[3];
sx q[3];
rz(-1.1464351) q[3];
sx q[3];
rz(-2.3800338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4789751) q[2];
sx q[2];
rz(-2.0765897) q[2];
sx q[2];
rz(0.62620658) q[2];
rz(-0.19691697) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(1.3885434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597252) q[0];
sx q[0];
rz(-1.4511061) q[0];
sx q[0];
rz(0.054280601) q[0];
rz(2.718603) q[1];
sx q[1];
rz(-1.3754247) q[1];
sx q[1];
rz(-1.8064226) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30935449) q[0];
sx q[0];
rz(-2.1119499) q[0];
sx q[0];
rz(-2.4914129) q[0];
rz(-0.35690709) q[2];
sx q[2];
rz(-2.1015328) q[2];
sx q[2];
rz(-2.4414947) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.74354913) q[1];
sx q[1];
rz(-0.7760007) q[1];
sx q[1];
rz(-2.3659124) q[1];
rz(3.0729592) q[3];
sx q[3];
rz(-1.4331927) q[3];
sx q[3];
rz(0.28561628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6173031) q[2];
sx q[2];
rz(-1.4681939) q[2];
sx q[2];
rz(2.7900901) q[2];
rz(-1.7678123) q[3];
sx q[3];
rz(-2.6280554) q[3];
sx q[3];
rz(2.8721749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3897301) q[0];
sx q[0];
rz(-2.8808012) q[0];
sx q[0];
rz(2.3824298) q[0];
rz(2.0579386) q[1];
sx q[1];
rz(-1.5307129) q[1];
sx q[1];
rz(-2.0738475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46684346) q[0];
sx q[0];
rz(-0.7848133) q[0];
sx q[0];
rz(0.87133566) q[0];
x q[1];
rz(0.2944417) q[2];
sx q[2];
rz(-1.8368372) q[2];
sx q[2];
rz(2.5584115) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8651004) q[1];
sx q[1];
rz(-2.4363764) q[1];
sx q[1];
rz(2.221285) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33371933) q[3];
sx q[3];
rz(-2.2401143) q[3];
sx q[3];
rz(-0.02726083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7119673) q[2];
sx q[2];
rz(-1.9316614) q[2];
sx q[2];
rz(0.21027002) q[2];
rz(-0.78768864) q[3];
sx q[3];
rz(-1.382788) q[3];
sx q[3];
rz(-2.6113966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44337153) q[0];
sx q[0];
rz(-2.5826695) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(-1.632985) q[1];
sx q[1];
rz(-1.2651545) q[1];
sx q[1];
rz(-1.9427585) q[1];
rz(-2.7349823) q[2];
sx q[2];
rz(-0.92401531) q[2];
sx q[2];
rz(0.32377908) q[2];
rz(1.8456712) q[3];
sx q[3];
rz(-1.7206921) q[3];
sx q[3];
rz(-0.10573798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
