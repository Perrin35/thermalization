OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(0.16790976) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(-2.8462703) q[1];
sx q[1];
rz(0.056161031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6457155) q[0];
sx q[0];
rz(-2.0464532) q[0];
sx q[0];
rz(0.23947421) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2916318) q[2];
sx q[2];
rz(-2.4764875) q[2];
sx q[2];
rz(-1.4814188) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4178884) q[1];
sx q[1];
rz(-2.1065518) q[1];
sx q[1];
rz(0.31236155) q[1];
x q[2];
rz(1.5264741) q[3];
sx q[3];
rz(-1.830415) q[3];
sx q[3];
rz(-2.0093902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7636259) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(1.1928605) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(-2.7584934) q[3];
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
rz(0.52779657) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(1.312785) q[0];
rz(-0.20547543) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(1.1516494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1412927) q[0];
sx q[0];
rz(-2.2465758) q[0];
sx q[0];
rz(-0.38210259) q[0];
rz(-pi) q[1];
rz(2.5580514) q[2];
sx q[2];
rz(-1.984664) q[2];
sx q[2];
rz(-1.5838503) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.68576751) q[1];
sx q[1];
rz(-2.0058504) q[1];
sx q[1];
rz(-2.0352092) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6395256) q[3];
sx q[3];
rz(-1.9126529) q[3];
sx q[3];
rz(-1.348192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1318704) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(-0.37718537) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(-2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31056988) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(-3.1047399) q[0];
rz(2.316078) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(-3.085014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6037613) q[0];
sx q[0];
rz(-2.1268401) q[0];
sx q[0];
rz(-2.6105196) q[0];
rz(0.25643202) q[2];
sx q[2];
rz(-1.7000546) q[2];
sx q[2];
rz(-1.749922) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.50443447) q[1];
sx q[1];
rz(-1.27099) q[1];
sx q[1];
rz(-1.7791041) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4967381) q[3];
sx q[3];
rz(-1.8861024) q[3];
sx q[3];
rz(0.23526084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8686707) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(2.2154714) q[2];
rz(0.55666322) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(-1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91519231) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(-0.74209374) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(2.6779968) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0915506) q[0];
sx q[0];
rz(-2.3766962) q[0];
sx q[0];
rz(-1.1372304) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.17642994) q[2];
sx q[2];
rz(-2.8039805) q[2];
sx q[2];
rz(-0.38844973) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1222893) q[1];
sx q[1];
rz(-2.3483854) q[1];
sx q[1];
rz(-1.8122458) q[1];
rz(-pi) q[2];
rz(0.47834088) q[3];
sx q[3];
rz(-1.0676749) q[3];
sx q[3];
rz(-2.2174045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3670369) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(3.0920933) q[2];
rz(-3.0130623) q[3];
sx q[3];
rz(-1.5934207) q[3];
sx q[3];
rz(0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-0.13609919) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(0.29770011) q[0];
rz(-0.4822576) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(0.94435) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56086841) q[0];
sx q[0];
rz(-1.7797911) q[0];
sx q[0];
rz(-3.0955549) q[0];
rz(-pi) q[1];
rz(-0.77563939) q[2];
sx q[2];
rz(-1.8619985) q[2];
sx q[2];
rz(-1.2736543) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9152865) q[1];
sx q[1];
rz(-1.5729135) q[1];
sx q[1];
rz(-1.6318984) q[1];
x q[2];
rz(0.51575757) q[3];
sx q[3];
rz(-0.48832794) q[3];
sx q[3];
rz(-2.8749136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8828316) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(-3.0333701) q[2];
rz(-0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.7047983) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(-0.58445245) q[0];
rz(-2.2553518) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(0.054919682) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9960105) q[0];
sx q[0];
rz(-1.3025563) q[0];
sx q[0];
rz(0.085573816) q[0];
rz(-pi) q[1];
rz(-1.6266277) q[2];
sx q[2];
rz(-0.32368127) q[2];
sx q[2];
rz(-1.2602381) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0262895) q[1];
sx q[1];
rz(-0.68740986) q[1];
sx q[1];
rz(2.5251212) q[1];
rz(-pi) q[2];
rz(2.218607) q[3];
sx q[3];
rz(-0.69159782) q[3];
sx q[3];
rz(1.637961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3871258) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(2.1248655) q[2];
rz(-0.54404849) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465437) q[0];
sx q[0];
rz(-0.98452079) q[0];
sx q[0];
rz(2.8570535) q[0];
rz(-2.1971205) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(-2.231266) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7401687) q[0];
sx q[0];
rz(-1.6356042) q[0];
sx q[0];
rz(3.0868953) q[0];
rz(-pi) q[1];
rz(-2.2482713) q[2];
sx q[2];
rz(-0.51572323) q[2];
sx q[2];
rz(2.233778) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0890761) q[1];
sx q[1];
rz(-2.0564582) q[1];
sx q[1];
rz(1.3658701) q[1];
rz(-pi) q[2];
rz(-2.5038239) q[3];
sx q[3];
rz(-2.31156) q[3];
sx q[3];
rz(-2.4142767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7515144) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(-2.2195623) q[2];
rz(2.5743124) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(-1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2475125) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(3.0122053) q[0];
rz(-2.5091876) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(0.30050373) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26586543) q[0];
sx q[0];
rz(-0.89571307) q[0];
sx q[0];
rz(0.48203326) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6408765) q[2];
sx q[2];
rz(-0.87101988) q[2];
sx q[2];
rz(-1.7388294) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0771675) q[1];
sx q[1];
rz(-2.7556813) q[1];
sx q[1];
rz(0.47972958) q[1];
x q[2];
rz(2.825533) q[3];
sx q[3];
rz(-2.3028767) q[3];
sx q[3];
rz(2.4162606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5552716) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(-0.78197455) q[2];
rz(-0.55109763) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5732116) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(-0.12776275) q[0];
rz(2.5993775) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(-2.382747) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7153873) q[0];
sx q[0];
rz(-1.1466768) q[0];
sx q[0];
rz(-0.018652648) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3779638) q[2];
sx q[2];
rz(-2.170993) q[2];
sx q[2];
rz(0.45229518) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85941852) q[1];
sx q[1];
rz(-2.4255883) q[1];
sx q[1];
rz(0.23846682) q[1];
rz(-pi) q[2];
x q[2];
rz(3.127029) q[3];
sx q[3];
rz(-2.2363538) q[3];
sx q[3];
rz(-1.927782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0163429) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(2.6515567) q[2];
rz(-1.7193433) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35995099) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(3.066257) q[0];
rz(2.244859) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(-2.5316701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88203726) q[0];
sx q[0];
rz(-1.8543188) q[0];
sx q[0];
rz(0.77772227) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8462734) q[2];
sx q[2];
rz(-1.2071949) q[2];
sx q[2];
rz(1.7737349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1249485) q[1];
sx q[1];
rz(-0.84608191) q[1];
sx q[1];
rz(-1.8925341) q[1];
rz(-pi) q[2];
rz(1.0697332) q[3];
sx q[3];
rz(-2.2483629) q[3];
sx q[3];
rz(2.5354405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23218368) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(-0.71371901) q[2];
rz(-0.37832007) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338035) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(-2.4907885) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(1.7667608) q[2];
sx q[2];
rz(-1.5407731) q[2];
sx q[2];
rz(-2.4827448) q[2];
rz(-1.1013423) q[3];
sx q[3];
rz(-2.6228842) q[3];
sx q[3];
rz(1.7026671) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];