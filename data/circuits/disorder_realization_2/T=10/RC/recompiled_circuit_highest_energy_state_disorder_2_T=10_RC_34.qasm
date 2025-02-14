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
rz(2.958137) q[0];
sx q[0];
rz(-2.3975211) q[0];
sx q[0];
rz(-2.0896572) q[0];
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
rz(-1.6010701) q[0];
sx q[0];
rz(-0.94694505) q[0];
sx q[0];
rz(-3.0649351) q[0];
rz(-2.6296205) q[2];
sx q[2];
rz(-2.3880771) q[2];
sx q[2];
rz(1.7640424) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4595653) q[1];
sx q[1];
rz(-0.50217705) q[1];
sx q[1];
rz(-2.9952496) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4450847) q[3];
sx q[3];
rz(-1.4918431) q[3];
sx q[3];
rz(-0.95206184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.28213349) q[2];
sx q[2];
rz(-0.30974516) q[2];
sx q[2];
rz(-2.0852883) q[2];
rz(3.1193962) q[3];
sx q[3];
rz(-2.3789417) q[3];
sx q[3];
rz(1.4200042) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9098772) q[0];
sx q[0];
rz(-0.48119369) q[0];
sx q[0];
rz(2.7114482) q[0];
rz(-3.0138956) q[1];
sx q[1];
rz(-1.1558497) q[1];
sx q[1];
rz(1.4375623) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6352954) q[0];
sx q[0];
rz(-0.40182913) q[0];
sx q[0];
rz(0.70962972) q[0];
rz(-2.2980437) q[2];
sx q[2];
rz(-3.0650716) q[2];
sx q[2];
rz(0.25329548) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.592652) q[1];
sx q[1];
rz(-0.90528622) q[1];
sx q[1];
rz(-1.1489465) q[1];
rz(-pi) q[2];
rz(2.6019179) q[3];
sx q[3];
rz(-1.7575348) q[3];
sx q[3];
rz(2.5591171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.11640707) q[2];
sx q[2];
rz(-1.4613287) q[2];
sx q[2];
rz(2.6961668) q[2];
rz(-2.7323501) q[3];
sx q[3];
rz(-1.0990812) q[3];
sx q[3];
rz(-0.87944952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55260783) q[0];
sx q[0];
rz(-0.092985066) q[0];
sx q[0];
rz(0.71480042) q[0];
rz(-1.0572761) q[1];
sx q[1];
rz(-2.7209268) q[1];
sx q[1];
rz(-0.32726273) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11571685) q[0];
sx q[0];
rz(-2.1109606) q[0];
sx q[0];
rz(2.9469423) q[0];
rz(0.94560854) q[2];
sx q[2];
rz(-2.473712) q[2];
sx q[2];
rz(-0.57648522) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.63673692) q[1];
sx q[1];
rz(-1.5914006) q[1];
sx q[1];
rz(1.7223174) q[1];
rz(-pi) q[2];
rz(2.8873575) q[3];
sx q[3];
rz(-1.2336858) q[3];
sx q[3];
rz(1.997662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0032234) q[2];
sx q[2];
rz(-2.1681163) q[2];
sx q[2];
rz(-1.6376015) q[2];
rz(1.2184294) q[3];
sx q[3];
rz(-2.7259493) q[3];
sx q[3];
rz(0.98201069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0686491) q[0];
sx q[0];
rz(-1.8953841) q[0];
sx q[0];
rz(-0.15596998) q[0];
rz(-2.7929557) q[1];
sx q[1];
rz(-2.5376384) q[1];
sx q[1];
rz(1.9085931) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6207244) q[0];
sx q[0];
rz(-1.8308795) q[0];
sx q[0];
rz(-3.0164032) q[0];
rz(-pi) q[1];
rz(-1.2945588) q[2];
sx q[2];
rz(-0.66670185) q[2];
sx q[2];
rz(1.897097) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.33267) q[1];
sx q[1];
rz(-0.70007174) q[1];
sx q[1];
rz(2.0168522) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[1];
rz(-1.4431346) q[2];
sx q[2];
rz(-0.036018697) q[2];
sx q[2];
rz(-1.5114463) q[2];
rz(-1.1394507) q[3];
sx q[3];
rz(-1.3316493) q[3];
sx q[3];
rz(-0.99036923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.776942) q[0];
sx q[0];
rz(-0.41780892) q[0];
sx q[0];
rz(-1.9885709) q[0];
rz(-0.54368377) q[1];
sx q[1];
rz(-1.8749571) q[1];
sx q[1];
rz(2.8797454) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42191045) q[0];
sx q[0];
rz(-1.6043264) q[0];
sx q[0];
rz(-1.6543341) q[0];
rz(0.33436454) q[2];
sx q[2];
rz(-1.0131256) q[2];
sx q[2];
rz(1.1209436) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10298577) q[1];
sx q[1];
rz(-2.1521882) q[1];
sx q[1];
rz(1.8379148) q[1];
rz(-pi) q[2];
rz(-0.20528593) q[3];
sx q[3];
rz(-2.3488099) q[3];
sx q[3];
rz(-2.1403811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5220945) q[2];
sx q[2];
rz(-2.5358989) q[2];
sx q[2];
rz(-1.4052793) q[2];
rz(2.3960579) q[3];
sx q[3];
rz(-1.2657974) q[3];
sx q[3];
rz(0.94329992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(3.140542) q[0];
sx q[0];
rz(-3.0241835) q[0];
sx q[0];
rz(-0.35181272) q[0];
rz(0.44031269) q[1];
sx q[1];
rz(-1.3654717) q[1];
sx q[1];
rz(0.17131677) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3765613) q[0];
sx q[0];
rz(-2.2269571) q[0];
sx q[0];
rz(1.2622467) q[0];
x q[1];
rz(2.3883853) q[2];
sx q[2];
rz(-1.79617) q[2];
sx q[2];
rz(-1.6373375) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79119392) q[1];
sx q[1];
rz(-1.5916675) q[1];
sx q[1];
rz(-2.4430165) q[1];
rz(-2.3790509) q[3];
sx q[3];
rz(-2.2665215) q[3];
sx q[3];
rz(-0.97555893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.064284023) q[2];
sx q[2];
rz(-1.372154) q[2];
sx q[2];
rz(2.1997931) q[2];
rz(-1.9866379) q[3];
sx q[3];
rz(-2.933511) q[3];
sx q[3];
rz(-1.8010767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6963541) q[0];
sx q[0];
rz(-2.0360763) q[0];
sx q[0];
rz(-0.0048986991) q[0];
rz(-0.44662961) q[1];
sx q[1];
rz(-0.59372562) q[1];
sx q[1];
rz(-0.68797025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5410091) q[0];
sx q[0];
rz(-2.2760411) q[0];
sx q[0];
rz(-1.7974822) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6737988) q[2];
sx q[2];
rz(-0.67343283) q[2];
sx q[2];
rz(1.3040257) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1157284) q[1];
sx q[1];
rz(-2.0021571) q[1];
sx q[1];
rz(-0.23484767) q[1];
rz(0.68617188) q[3];
sx q[3];
rz(-0.99643512) q[3];
sx q[3];
rz(-2.2597093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.530431) q[2];
sx q[2];
rz(-2.7335584) q[2];
sx q[2];
rz(-0.92998695) q[2];
rz(-1.9200578) q[3];
sx q[3];
rz(-1.113021) q[3];
sx q[3];
rz(-0.59337029) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68194836) q[0];
sx q[0];
rz(-1.3197897) q[0];
sx q[0];
rz(-3.0514858) q[0];
rz(1.7168761) q[1];
sx q[1];
rz(-0.63980353) q[1];
sx q[1];
rz(2.7065014) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85283576) q[0];
sx q[0];
rz(-2.618874) q[0];
sx q[0];
rz(-2.8466892) q[0];
rz(-pi) q[1];
rz(-0.22771671) q[2];
sx q[2];
rz(-0.68074742) q[2];
sx q[2];
rz(-1.6344223) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.69999483) q[1];
sx q[1];
rz(-0.95788237) q[1];
sx q[1];
rz(0.15721486) q[1];
x q[2];
rz(-3.1396542) q[3];
sx q[3];
rz(-0.42436212) q[3];
sx q[3];
rz(-0.75943179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4789751) q[2];
sx q[2];
rz(-2.0765897) q[2];
sx q[2];
rz(-0.62620658) q[2];
rz(-2.9446757) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(-1.3885434) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58186746) q[0];
sx q[0];
rz(-1.6904866) q[0];
sx q[0];
rz(0.054280601) q[0];
rz(2.718603) q[1];
sx q[1];
rz(-1.7661679) q[1];
sx q[1];
rz(1.8064226) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88807073) q[0];
sx q[0];
rz(-1.0253064) q[0];
sx q[0];
rz(-2.2175199) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7846856) q[2];
sx q[2];
rz(-1.0400598) q[2];
sx q[2];
rz(-0.70009795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21653342) q[1];
sx q[1];
rz(-1.0581985) q[1];
sx q[1];
rz(2.53043) q[1];
rz(1.7087206) q[3];
sx q[3];
rz(-1.63878) q[3];
sx q[3];
rz(1.2757511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6173031) q[2];
sx q[2];
rz(-1.6733988) q[2];
sx q[2];
rz(-2.7900901) q[2];
rz(-1.3737804) q[3];
sx q[3];
rz(-2.6280554) q[3];
sx q[3];
rz(0.26941776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7518625) q[0];
sx q[0];
rz(-2.8808012) q[0];
sx q[0];
rz(0.75916284) q[0];
rz(-1.0836541) q[1];
sx q[1];
rz(-1.5307129) q[1];
sx q[1];
rz(-2.0738475) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7368374) q[0];
sx q[0];
rz(-2.1421297) q[0];
sx q[0];
rz(-2.570117) q[0];
x q[1];
rz(-2.3876486) q[2];
sx q[2];
rz(-2.747376) q[2];
sx q[2];
rz(-2.8682402) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8651004) q[1];
sx q[1];
rz(-2.4363764) q[1];
sx q[1];
rz(-0.92030763) q[1];
rz(0.87369793) q[3];
sx q[3];
rz(-1.3109968) q[3];
sx q[3];
rz(1.3316621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4296253) q[2];
sx q[2];
rz(-1.9316614) q[2];
sx q[2];
rz(-0.21027002) q[2];
rz(2.353904) q[3];
sx q[3];
rz(-1.7588047) q[3];
sx q[3];
rz(2.6113966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44337153) q[0];
sx q[0];
rz(-0.5589232) q[0];
sx q[0];
rz(1.9236175) q[0];
rz(-1.5086077) q[1];
sx q[1];
rz(-1.8764381) q[1];
sx q[1];
rz(1.1988342) q[1];
rz(2.2589113) q[2];
sx q[2];
rz(-1.8919049) q[2];
sx q[2];
rz(1.6406825) q[2];
rz(-0.15565025) q[3];
sx q[3];
rz(-1.8425103) q[3];
sx q[3];
rz(1.507148) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
