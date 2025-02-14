OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2321229) q[0];
sx q[0];
rz(-0.98804086) q[0];
sx q[0];
rz(-2.7890132) q[0];
rz(-0.70631385) q[1];
sx q[1];
rz(-0.7847473) q[1];
sx q[1];
rz(0.027332505) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77307149) q[0];
sx q[0];
rz(-1.3734682) q[0];
sx q[0];
rz(2.7993566) q[0];
x q[1];
rz(2.4800221) q[2];
sx q[2];
rz(-2.9230766) q[2];
sx q[2];
rz(-2.0486346) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.24522745) q[1];
sx q[1];
rz(-1.7385529) q[1];
sx q[1];
rz(1.5390525) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23689278) q[3];
sx q[3];
rz(-1.4290571) q[3];
sx q[3];
rz(0.6976034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0722384) q[2];
sx q[2];
rz(-1.8202929) q[2];
sx q[2];
rz(1.5864774) q[2];
rz(2.41411) q[3];
sx q[3];
rz(-2.1860217) q[3];
sx q[3];
rz(2.6684842) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.469406) q[0];
sx q[0];
rz(-1.3453901) q[0];
sx q[0];
rz(2.192002) q[0];
rz(2.8166215) q[1];
sx q[1];
rz(-1.800671) q[1];
sx q[1];
rz(-0.51345888) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0605704) q[0];
sx q[0];
rz(-2.1784884) q[0];
sx q[0];
rz(-0.49852235) q[0];
rz(-pi) q[1];
rz(1.9302692) q[2];
sx q[2];
rz(-2.3811334) q[2];
sx q[2];
rz(2.9414842) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.70043515) q[1];
sx q[1];
rz(-1.0714021) q[1];
sx q[1];
rz(2.5446289) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6504032) q[3];
sx q[3];
rz(-2.0955387) q[3];
sx q[3];
rz(0.82651807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.010633858) q[2];
sx q[2];
rz(-1.3291239) q[2];
sx q[2];
rz(1.0884253) q[2];
rz(-2.4736577) q[3];
sx q[3];
rz(-0.1440983) q[3];
sx q[3];
rz(1.7264504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82069355) q[0];
sx q[0];
rz(-2.2405393) q[0];
sx q[0];
rz(-1.2574842) q[0];
rz(0.084818689) q[1];
sx q[1];
rz(-0.091400472) q[1];
sx q[1];
rz(0.65748293) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5710058) q[0];
sx q[0];
rz(-2.0181542) q[0];
sx q[0];
rz(-1.7156832) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0107747) q[2];
sx q[2];
rz(-2.0891671) q[2];
sx q[2];
rz(0.033998297) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7967826) q[1];
sx q[1];
rz(-1.0940897) q[1];
sx q[1];
rz(1.7631129) q[1];
rz(-pi) q[2];
rz(1.182193) q[3];
sx q[3];
rz(-2.4425584) q[3];
sx q[3];
rz(0.70632589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76454863) q[2];
sx q[2];
rz(-0.42082861) q[2];
sx q[2];
rz(1.8910889) q[2];
rz(-1.4416384) q[3];
sx q[3];
rz(-1.3911894) q[3];
sx q[3];
rz(-0.82956782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4568951) q[0];
sx q[0];
rz(-0.95624113) q[0];
sx q[0];
rz(-0.60234219) q[0];
rz(-2.1777228) q[1];
sx q[1];
rz(-0.78097051) q[1];
sx q[1];
rz(-0.22183713) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6312799) q[0];
sx q[0];
rz(-1.5595734) q[0];
sx q[0];
rz(1.5785286) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1603951) q[2];
sx q[2];
rz(-2.5735435) q[2];
sx q[2];
rz(-2.8063453) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.027699359) q[1];
sx q[1];
rz(-2.4287819) q[1];
sx q[1];
rz(-0.63473397) q[1];
rz(-pi) q[2];
rz(2.6915914) q[3];
sx q[3];
rz(-1.5061989) q[3];
sx q[3];
rz(-0.76417506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43231371) q[2];
sx q[2];
rz(-1.1608492) q[2];
sx q[2];
rz(-0.088689001) q[2];
rz(0.49897075) q[3];
sx q[3];
rz(-1.5191398) q[3];
sx q[3];
rz(-0.078908198) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13224193) q[0];
sx q[0];
rz(-3.0920588) q[0];
sx q[0];
rz(2.2392739) q[0];
rz(-2.8514476) q[1];
sx q[1];
rz(-0.90466181) q[1];
sx q[1];
rz(-1.9754999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9525495) q[0];
sx q[0];
rz(-1.4092085) q[0];
sx q[0];
rz(2.5717006) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3275215) q[2];
sx q[2];
rz(-1.8406879) q[2];
sx q[2];
rz(-0.44651595) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.13671) q[1];
sx q[1];
rz(-0.45082475) q[1];
sx q[1];
rz(-0.87597998) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5447292) q[3];
sx q[3];
rz(-1.8223624) q[3];
sx q[3];
rz(0.1840011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0151998) q[2];
sx q[2];
rz(-0.79605278) q[2];
sx q[2];
rz(-1.2882721) q[2];
rz(-1.0334233) q[3];
sx q[3];
rz(-1.1181592) q[3];
sx q[3];
rz(-1.6667295) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.425659) q[0];
sx q[0];
rz(-0.10361828) q[0];
sx q[0];
rz(-0.15956751) q[0];
rz(-0.60699925) q[1];
sx q[1];
rz(-1.9885352) q[1];
sx q[1];
rz(-2.0807696) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.589405) q[0];
sx q[0];
rz(-0.97597835) q[0];
sx q[0];
rz(2.908559) q[0];
rz(-0.26788864) q[2];
sx q[2];
rz(-1.3527591) q[2];
sx q[2];
rz(1.5341369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6276803) q[1];
sx q[1];
rz(-2.5315171) q[1];
sx q[1];
rz(1.3327636) q[1];
x q[2];
rz(0.14935519) q[3];
sx q[3];
rz(-0.58404446) q[3];
sx q[3];
rz(-0.67549878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1242421) q[2];
sx q[2];
rz(-0.81581798) q[2];
sx q[2];
rz(0.38786495) q[2];
rz(-1.806949) q[3];
sx q[3];
rz(-2.4692061) q[3];
sx q[3];
rz(2.1147125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89011985) q[0];
sx q[0];
rz(-2.1595182) q[0];
sx q[0];
rz(0.89333308) q[0];
rz(-2.6003301) q[1];
sx q[1];
rz(-0.90507871) q[1];
sx q[1];
rz(2.0792686) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5389299) q[0];
sx q[0];
rz(-1.8315803) q[0];
sx q[0];
rz(-2.3988924) q[0];
rz(-pi) q[1];
rz(-0.063060305) q[2];
sx q[2];
rz(-0.82411924) q[2];
sx q[2];
rz(-2.7450739) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4121288) q[1];
sx q[1];
rz(-2.1264014) q[1];
sx q[1];
rz(1.5430862) q[1];
rz(-2.9422308) q[3];
sx q[3];
rz(-2.3219675) q[3];
sx q[3];
rz(0.79341753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5541151) q[2];
sx q[2];
rz(-0.97439659) q[2];
sx q[2];
rz(-1.1678196) q[2];
rz(0.81985146) q[3];
sx q[3];
rz(-2.3746115) q[3];
sx q[3];
rz(0.64275536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5646566) q[0];
sx q[0];
rz(-2.7695203) q[0];
sx q[0];
rz(1.8336953) q[0];
rz(1.4564184) q[1];
sx q[1];
rz(-1.5142199) q[1];
sx q[1];
rz(0.13052043) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3878138) q[0];
sx q[0];
rz(-0.31720933) q[0];
sx q[0];
rz(1.8075726) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0840262) q[2];
sx q[2];
rz(-2.3064838) q[2];
sx q[2];
rz(1.4722919) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.91192818) q[1];
sx q[1];
rz(-2.4897631) q[1];
sx q[1];
rz(-0.064464557) q[1];
rz(2.42584) q[3];
sx q[3];
rz(-2.0945315) q[3];
sx q[3];
rz(0.88818404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.1859583) q[2];
sx q[2];
rz(-0.83883494) q[2];
sx q[2];
rz(-1.4107417) q[2];
rz(-3.020982) q[3];
sx q[3];
rz(-1.6552304) q[3];
sx q[3];
rz(-1.6284778) q[3];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8677583) q[0];
sx q[0];
rz(-0.61536106) q[0];
sx q[0];
rz(-0.14061418) q[0];
rz(0.73453844) q[1];
sx q[1];
rz(-1.7081552) q[1];
sx q[1];
rz(-2.3908884) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.984042) q[0];
sx q[0];
rz(-0.37985248) q[0];
sx q[0];
rz(-0.67870875) q[0];
rz(-1.2241988) q[2];
sx q[2];
rz(-1.7945514) q[2];
sx q[2];
rz(2.4181197) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0486601) q[1];
sx q[1];
rz(-1.5092998) q[1];
sx q[1];
rz(-1.3912806) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4846913) q[3];
sx q[3];
rz(-0.70843452) q[3];
sx q[3];
rz(-1.3902067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1371896) q[2];
sx q[2];
rz(-0.69978324) q[2];
sx q[2];
rz(-2.3889551) q[2];
rz(-0.33958069) q[3];
sx q[3];
rz(-1.7759674) q[3];
sx q[3];
rz(0.021421758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.705377) q[0];
sx q[0];
rz(-0.74420539) q[0];
sx q[0];
rz(-0.63802737) q[0];
rz(-0.72884196) q[1];
sx q[1];
rz(-1.3325656) q[1];
sx q[1];
rz(-2.3400838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8014098) q[0];
sx q[0];
rz(-1.3687642) q[0];
sx q[0];
rz(2.3372746) q[0];
rz(-pi) q[1];
rz(-0.83790581) q[2];
sx q[2];
rz(-1.3896754) q[2];
sx q[2];
rz(-1.52219) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2406225) q[1];
sx q[1];
rz(-0.83604807) q[1];
sx q[1];
rz(-0.052608629) q[1];
rz(-2.2156694) q[3];
sx q[3];
rz(-1.7002673) q[3];
sx q[3];
rz(2.7162939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.30404299) q[2];
sx q[2];
rz(-1.3648698) q[2];
sx q[2];
rz(2.020828) q[2];
rz(2.4325727) q[3];
sx q[3];
rz(-2.6778383) q[3];
sx q[3];
rz(-1.2254865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7365702) q[0];
sx q[0];
rz(-2.6800192) q[0];
sx q[0];
rz(-0.17929684) q[0];
rz(-2.2364521) q[1];
sx q[1];
rz(-2.1642579) q[1];
sx q[1];
rz(-0.49107818) q[1];
rz(0.35781833) q[2];
sx q[2];
rz(-2.2651548) q[2];
sx q[2];
rz(1.6969778) q[2];
rz(-2.0427451) q[3];
sx q[3];
rz(-0.57335486) q[3];
sx q[3];
rz(-2.2344979) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
