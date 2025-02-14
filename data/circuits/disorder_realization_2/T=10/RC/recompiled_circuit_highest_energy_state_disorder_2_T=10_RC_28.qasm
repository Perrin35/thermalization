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
rz(0.19368859) q[1];
sx q[1];
rz(-0.63313484) q[1];
sx q[1];
rz(-2.5685891) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0664803) q[0];
sx q[0];
rz(-1.5085992) q[0];
sx q[0];
rz(0.94554995) q[0];
x q[1];
rz(-0.68555514) q[2];
sx q[2];
rz(-1.9126045) q[2];
sx q[2];
rz(2.5594001) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.51542789) q[1];
sx q[1];
rz(-2.0671131) q[1];
sx q[1];
rz(-1.650701) q[1];
x q[2];
rz(1.4680193) q[3];
sx q[3];
rz(-2.2647018) q[3];
sx q[3];
rz(0.68460195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.28213349) q[2];
sx q[2];
rz(-2.8318475) q[2];
sx q[2];
rz(1.0563043) q[2];
rz(3.1193962) q[3];
sx q[3];
rz(-0.76265097) q[3];
sx q[3];
rz(-1.4200042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9098772) q[0];
sx q[0];
rz(-2.660399) q[0];
sx q[0];
rz(-2.7114482) q[0];
rz(3.0138956) q[1];
sx q[1];
rz(-1.1558497) q[1];
sx q[1];
rz(1.7040303) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50629726) q[0];
sx q[0];
rz(-2.7397635) q[0];
sx q[0];
rz(2.4319629) q[0];
rz(-pi) q[1];
rz(2.2980437) q[2];
sx q[2];
rz(-0.076521046) q[2];
sx q[2];
rz(-2.8882972) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3900657) q[1];
sx q[1];
rz(-1.8987149) q[1];
sx q[1];
rz(-0.71050127) q[1];
rz(-2.7892465) q[3];
sx q[3];
rz(-0.56803136) q[3];
sx q[3];
rz(1.2888792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11640707) q[2];
sx q[2];
rz(-1.4613287) q[2];
sx q[2];
rz(-0.44542584) q[2];
rz(-2.7323501) q[3];
sx q[3];
rz(-1.0990812) q[3];
sx q[3];
rz(2.2621431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55260783) q[0];
sx q[0];
rz(-3.0486076) q[0];
sx q[0];
rz(2.4267922) q[0];
rz(-2.0843166) q[1];
sx q[1];
rz(-2.7209268) q[1];
sx q[1];
rz(0.32726273) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0258758) q[0];
sx q[0];
rz(-2.1109606) q[0];
sx q[0];
rz(-0.19465036) q[0];
rz(2.1398323) q[2];
sx q[2];
rz(-1.941701) q[2];
sx q[2];
rz(-0.47874622) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93720528) q[1];
sx q[1];
rz(-1.722285) q[1];
sx q[1];
rz(3.1207496) q[1];
x q[2];
rz(-0.9483665) q[3];
sx q[3];
rz(-2.7223058) q[3];
sx q[3];
rz(0.47800999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1383692) q[2];
sx q[2];
rz(-0.97347632) q[2];
sx q[2];
rz(1.5039911) q[2];
rz(1.9231632) q[3];
sx q[3];
rz(-0.4156433) q[3];
sx q[3];
rz(-2.159582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0686491) q[0];
sx q[0];
rz(-1.8953841) q[0];
sx q[0];
rz(0.15596998) q[0];
rz(-2.7929557) q[1];
sx q[1];
rz(-0.60395423) q[1];
sx q[1];
rz(1.2329996) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1656146) q[0];
sx q[0];
rz(-2.8535643) q[0];
sx q[0];
rz(-2.0095129) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2945588) q[2];
sx q[2];
rz(-2.4748908) q[2];
sx q[2];
rz(-1.897097) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2501331) q[1];
sx q[1];
rz(-2.1910408) q[1];
sx q[1];
rz(-0.34858443) q[1];
rz(1.5855012) q[3];
sx q[3];
rz(-1.8222408) q[3];
sx q[3];
rz(-0.78212839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6984581) q[2];
sx q[2];
rz(-0.036018697) q[2];
sx q[2];
rz(1.6301463) q[2];
rz(2.002142) q[3];
sx q[3];
rz(-1.3316493) q[3];
sx q[3];
rz(2.1512234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.776942) q[0];
sx q[0];
rz(-2.7237837) q[0];
sx q[0];
rz(-1.9885709) q[0];
rz(-0.54368377) q[1];
sx q[1];
rz(-1.8749571) q[1];
sx q[1];
rz(2.8797454) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76806289) q[0];
sx q[0];
rz(-3.0515915) q[0];
sx q[0];
rz(1.1885719) q[0];
x q[1];
rz(-1.0864429) q[2];
sx q[2];
rz(-0.64099714) q[2];
sx q[2];
rz(-0.5400368) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56526977) q[1];
sx q[1];
rz(-2.5082631) q[1];
sx q[1];
rz(-0.38193462) q[1];
x q[2];
rz(2.35942) q[3];
sx q[3];
rz(-1.425079) q[3];
sx q[3];
rz(-2.4268933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5220945) q[2];
sx q[2];
rz(-0.60569373) q[2];
sx q[2];
rz(-1.4052793) q[2];
rz(-2.3960579) q[3];
sx q[3];
rz(-1.8757952) q[3];
sx q[3];
rz(-2.1982927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-2.7897799) q[0];
rz(-2.70128) q[1];
sx q[1];
rz(-1.7761209) q[1];
sx q[1];
rz(-0.17131677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2835944) q[0];
sx q[0];
rz(-2.4263315) q[0];
sx q[0];
rz(-2.7659228) q[0];
rz(-1.8753042) q[2];
sx q[2];
rz(-0.84104462) q[2];
sx q[2];
rz(2.8685121) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80444634) q[1];
sx q[1];
rz(-2.4427572) q[1];
sx q[1];
rz(0.032445907) q[1];
rz(-pi) q[2];
rz(-0.71368696) q[3];
sx q[3];
rz(-1.011935) q[3];
sx q[3];
rz(-1.9969459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.064284023) q[2];
sx q[2];
rz(-1.372154) q[2];
sx q[2];
rz(-0.94179955) q[2];
rz(-1.9866379) q[3];
sx q[3];
rz(-2.933511) q[3];
sx q[3];
rz(-1.8010767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.6963541) q[0];
sx q[0];
rz(-2.0360763) q[0];
sx q[0];
rz(3.136694) q[0];
rz(2.694963) q[1];
sx q[1];
rz(-0.59372562) q[1];
sx q[1];
rz(-0.68797025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82179994) q[0];
sx q[0];
rz(-1.3988136) q[0];
sx q[0];
rz(-0.71806192) q[0];
x q[1];
rz(2.5227658) q[2];
sx q[2];
rz(-1.2857253) q[2];
sx q[2];
rz(0.10933354) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5060978) q[1];
sx q[1];
rz(-2.6539972) q[1];
sx q[1];
rz(-1.1027085) q[1];
x q[2];
rz(0.68617188) q[3];
sx q[3];
rz(-2.1451575) q[3];
sx q[3];
rz(2.2597093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61116162) q[2];
sx q[2];
rz(-2.7335584) q[2];
sx q[2];
rz(-0.92998695) q[2];
rz(1.9200578) q[3];
sx q[3];
rz(-1.113021) q[3];
sx q[3];
rz(-2.5482224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4596443) q[0];
sx q[0];
rz(-1.3197897) q[0];
sx q[0];
rz(0.090106877) q[0];
rz(1.7168761) q[1];
sx q[1];
rz(-2.5017891) q[1];
sx q[1];
rz(0.43509126) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85283576) q[0];
sx q[0];
rz(-2.618874) q[0];
sx q[0];
rz(-0.29490348) q[0];
rz(-0.66799156) q[2];
sx q[2];
rz(-1.7133623) q[2];
sx q[2];
rz(-0.1145471) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96174091) q[1];
sx q[1];
rz(-1.6992178) q[1];
sx q[1];
rz(-0.95203103) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7172312) q[3];
sx q[3];
rz(-1.5715944) q[3];
sx q[3];
rz(-0.80959807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4789751) q[2];
sx q[2];
rz(-2.0765897) q[2];
sx q[2];
rz(2.5153861) q[2];
rz(2.9446757) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(-1.7530493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58186746) q[0];
sx q[0];
rz(-1.6904866) q[0];
sx q[0];
rz(-0.054280601) q[0];
rz(2.718603) q[1];
sx q[1];
rz(-1.7661679) q[1];
sx q[1];
rz(-1.3351701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88807073) q[0];
sx q[0];
rz(-1.0253064) q[0];
sx q[0];
rz(2.2175199) q[0];
rz(-1.0112052) q[2];
sx q[2];
rz(-1.8768684) q[2];
sx q[2];
rz(-1.0572421) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9250592) q[1];
sx q[1];
rz(-1.0581985) q[1];
sx q[1];
rz(2.53043) q[1];
x q[2];
rz(3.0729592) q[3];
sx q[3];
rz(-1.4331927) q[3];
sx q[3];
rz(0.28561628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52428952) q[2];
sx q[2];
rz(-1.4681939) q[2];
sx q[2];
rz(0.35150251) q[2];
rz(1.3737804) q[3];
sx q[3];
rz(-0.51353729) q[3];
sx q[3];
rz(-2.8721749) q[3];
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
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7518625) q[0];
sx q[0];
rz(-0.26079145) q[0];
sx q[0];
rz(-0.75916284) q[0];
rz(2.0579386) q[1];
sx q[1];
rz(-1.5307129) q[1];
sx q[1];
rz(1.0677451) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6747492) q[0];
sx q[0];
rz(-2.3567794) q[0];
sx q[0];
rz(2.270257) q[0];
x q[1];
rz(-1.8482089) q[2];
sx q[2];
rz(-1.8545863) q[2];
sx q[2];
rz(-2.0744155) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.23087654) q[1];
sx q[1];
rz(-1.1674122) q[1];
sx q[1];
rz(-2.1662) q[1];
x q[2];
rz(-2.2678947) q[3];
sx q[3];
rz(-1.3109968) q[3];
sx q[3];
rz(-1.8099305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4296253) q[2];
sx q[2];
rz(-1.2099313) q[2];
sx q[2];
rz(-0.21027002) q[2];
rz(2.353904) q[3];
sx q[3];
rz(-1.382788) q[3];
sx q[3];
rz(0.53019607) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(2.0532578) q[2];
sx q[2];
rz(-2.3934622) q[2];
sx q[2];
rz(-0.29665034) q[2];
rz(0.15565025) q[3];
sx q[3];
rz(-1.2990824) q[3];
sx q[3];
rz(-1.6344447) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
