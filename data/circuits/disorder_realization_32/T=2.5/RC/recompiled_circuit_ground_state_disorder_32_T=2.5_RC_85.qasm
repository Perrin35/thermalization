OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6736421) q[0];
sx q[0];
rz(-1.8704432) q[0];
sx q[0];
rz(-0.76572642) q[0];
rz(2.7453121) q[1];
sx q[1];
rz(-3.0488465) q[1];
sx q[1];
rz(-0.5527817) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3512289) q[0];
sx q[0];
rz(-1.1811331) q[0];
sx q[0];
rz(-2.0084951) q[0];
rz(-pi) q[1];
rz(-2.2596006) q[2];
sx q[2];
rz(-2.246736) q[2];
sx q[2];
rz(1.9641528) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0243135) q[1];
sx q[1];
rz(-1.8872042) q[1];
sx q[1];
rz(0.98543075) q[1];
rz(-pi) q[2];
rz(2.1215274) q[3];
sx q[3];
rz(-2.4006776) q[3];
sx q[3];
rz(-1.7307841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7736241) q[2];
sx q[2];
rz(-1.5017193) q[2];
sx q[2];
rz(3.0527414) q[2];
rz(2.7446274) q[3];
sx q[3];
rz(-1.0395972) q[3];
sx q[3];
rz(-1.6840434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48001584) q[0];
sx q[0];
rz(-0.50742298) q[0];
sx q[0];
rz(-0.85451025) q[0];
rz(-2.5201733) q[1];
sx q[1];
rz(-2.8050551) q[1];
sx q[1];
rz(-2.0889166) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39256477) q[0];
sx q[0];
rz(-1.6727007) q[0];
sx q[0];
rz(0.13691957) q[0];
rz(2.3143155) q[2];
sx q[2];
rz(-0.57663871) q[2];
sx q[2];
rz(1.0232384) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4667644) q[1];
sx q[1];
rz(-0.75777868) q[1];
sx q[1];
rz(1.8653052) q[1];
rz(1.5731811) q[3];
sx q[3];
rz(-2.5376476) q[3];
sx q[3];
rz(-1.9238453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.66775995) q[2];
sx q[2];
rz(-0.86060539) q[2];
sx q[2];
rz(2.1916126) q[2];
rz(1.0085603) q[3];
sx q[3];
rz(-0.63298321) q[3];
sx q[3];
rz(2.8620201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94473332) q[0];
sx q[0];
rz(-1.9254528) q[0];
sx q[0];
rz(-1.7945633) q[0];
rz(1.0579717) q[1];
sx q[1];
rz(-0.24669138) q[1];
sx q[1];
rz(-0.11014858) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0824688) q[0];
sx q[0];
rz(-1.1131071) q[0];
sx q[0];
rz(2.7864866) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0588689) q[2];
sx q[2];
rz(-0.16943585) q[2];
sx q[2];
rz(-1.8878586) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5837142) q[1];
sx q[1];
rz(-1.4643351) q[1];
sx q[1];
rz(-2.540968) q[1];
rz(0.44656698) q[3];
sx q[3];
rz(-2.0166698) q[3];
sx q[3];
rz(1.879541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.18016732) q[2];
sx q[2];
rz(-1.578293) q[2];
sx q[2];
rz(3.0999198) q[2];
rz(-2.9336119) q[3];
sx q[3];
rz(-0.22170034) q[3];
sx q[3];
rz(2.8878133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5792849) q[0];
sx q[0];
rz(-1.8203745) q[0];
sx q[0];
rz(2.1606309) q[0];
rz(-2.7085069) q[1];
sx q[1];
rz(-1.667645) q[1];
sx q[1];
rz(1.6281737) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5034803) q[0];
sx q[0];
rz(-0.44027087) q[0];
sx q[0];
rz(-0.43285893) q[0];
x q[1];
rz(-1.5536521) q[2];
sx q[2];
rz(-2.4701256) q[2];
sx q[2];
rz(1.4070828) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.42964307) q[1];
sx q[1];
rz(-2.6338661) q[1];
sx q[1];
rz(-1.7832548) q[1];
rz(0.93928503) q[3];
sx q[3];
rz(-1.8016677) q[3];
sx q[3];
rz(0.38294562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5341586) q[2];
sx q[2];
rz(-0.51924339) q[2];
sx q[2];
rz(-2.329211) q[2];
rz(-1.2477929) q[3];
sx q[3];
rz(-1.2500074) q[3];
sx q[3];
rz(1.2468503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.3355376) q[0];
sx q[0];
rz(-2.1195109) q[0];
sx q[0];
rz(-1.442765) q[0];
rz(2.6834148) q[1];
sx q[1];
rz(-1.5135601) q[1];
sx q[1];
rz(-2.3853669) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6578015) q[0];
sx q[0];
rz(-1.4157802) q[0];
sx q[0];
rz(2.20138) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5014072) q[2];
sx q[2];
rz(-2.019503) q[2];
sx q[2];
rz(1.9404836) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9971265) q[1];
sx q[1];
rz(-1.4835852) q[1];
sx q[1];
rz(-3.0060796) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9147768) q[3];
sx q[3];
rz(-1.6046673) q[3];
sx q[3];
rz(-0.86465166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4592287) q[2];
sx q[2];
rz(-1.3882779) q[2];
sx q[2];
rz(0.96561042) q[2];
rz(-2.5718578) q[3];
sx q[3];
rz(-1.1252879) q[3];
sx q[3];
rz(1.4558571) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7392893) q[0];
sx q[0];
rz(-1.9434513) q[0];
sx q[0];
rz(-1.7560316) q[0];
rz(2.3784474) q[1];
sx q[1];
rz(-1.4475854) q[1];
sx q[1];
rz(0.10173434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83086936) q[0];
sx q[0];
rz(-0.95267696) q[0];
sx q[0];
rz(-0.11316664) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5237807) q[2];
sx q[2];
rz(-1.1883548) q[2];
sx q[2];
rz(-2.2161992) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6468874) q[1];
sx q[1];
rz(-0.27783074) q[1];
sx q[1];
rz(-1.7956514) q[1];
x q[2];
rz(0.88512986) q[3];
sx q[3];
rz(-2.1477774) q[3];
sx q[3];
rz(-2.5654405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10852854) q[2];
sx q[2];
rz(-1.6988924) q[2];
sx q[2];
rz(-0.70890439) q[2];
rz(2.3809643) q[3];
sx q[3];
rz(-1.9899188) q[3];
sx q[3];
rz(0.93543783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0660504) q[0];
sx q[0];
rz(-1.1533371) q[0];
sx q[0];
rz(0.72108889) q[0];
rz(0.9779633) q[1];
sx q[1];
rz(-2.5826192) q[1];
sx q[1];
rz(1.5616547) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52873351) q[0];
sx q[0];
rz(-1.4943143) q[0];
sx q[0];
rz(-1.6441109) q[0];
rz(-pi) q[1];
rz(-2.8992813) q[2];
sx q[2];
rz(-2.152395) q[2];
sx q[2];
rz(2.966449) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9813919) q[1];
sx q[1];
rz(-1.0202279) q[1];
sx q[1];
rz(-1.1070613) q[1];
rz(2.0098814) q[3];
sx q[3];
rz(-2.4222322) q[3];
sx q[3];
rz(0.92911014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28479031) q[2];
sx q[2];
rz(-0.52758104) q[2];
sx q[2];
rz(-1.6214726) q[2];
rz(2.704845) q[3];
sx q[3];
rz(-0.89478651) q[3];
sx q[3];
rz(0.53957087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.082551) q[0];
sx q[0];
rz(-0.83913791) q[0];
sx q[0];
rz(-2.3178597) q[0];
rz(-2.0027022) q[1];
sx q[1];
rz(-1.7887807) q[1];
sx q[1];
rz(-1.9653856) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19768077) q[0];
sx q[0];
rz(-1.0689387) q[0];
sx q[0];
rz(-3.0905064) q[0];
rz(-pi) q[1];
rz(-1.9763821) q[2];
sx q[2];
rz(-2.1422804) q[2];
sx q[2];
rz(1.6125082) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9261647) q[1];
sx q[1];
rz(-1.3326982) q[1];
sx q[1];
rz(1.560731) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.05925624) q[3];
sx q[3];
rz(-2.5496581) q[3];
sx q[3];
rz(-1.8979285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99522432) q[2];
sx q[2];
rz(-2.5547042) q[2];
sx q[2];
rz(0.5874908) q[2];
rz(2.7557709) q[3];
sx q[3];
rz(-0.84857517) q[3];
sx q[3];
rz(-2.2390656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0226444) q[0];
sx q[0];
rz(-0.97624874) q[0];
sx q[0];
rz(-0.60978419) q[0];
rz(-1.7976286) q[1];
sx q[1];
rz(-1.983843) q[1];
sx q[1];
rz(-2.9885898) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5822844) q[0];
sx q[0];
rz(-2.0366208) q[0];
sx q[0];
rz(1.5488173) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79522466) q[2];
sx q[2];
rz(-1.5780026) q[2];
sx q[2];
rz(1.1178655) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3045522) q[1];
sx q[1];
rz(-0.64230949) q[1];
sx q[1];
rz(1.6019956) q[1];
rz(-pi) q[2];
rz(2.6681247) q[3];
sx q[3];
rz(-2.2064379) q[3];
sx q[3];
rz(-1.1421695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0387705) q[2];
sx q[2];
rz(-1.5800579) q[2];
sx q[2];
rz(2.5193396) q[2];
rz(1.2004987) q[3];
sx q[3];
rz(-1.7222722) q[3];
sx q[3];
rz(2.5672369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822561) q[0];
sx q[0];
rz(-0.065956235) q[0];
sx q[0];
rz(-0.35520735) q[0];
rz(2.0390873) q[1];
sx q[1];
rz(-3.1246154) q[1];
sx q[1];
rz(0.65974081) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0416464) q[0];
sx q[0];
rz(-0.92586854) q[0];
sx q[0];
rz(1.7862595) q[0];
rz(0.12730403) q[2];
sx q[2];
rz(-1.9812756) q[2];
sx q[2];
rz(-2.8212027) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8568032) q[1];
sx q[1];
rz(-1.2592013) q[1];
sx q[1];
rz(1.9143399) q[1];
x q[2];
rz(-1.6361314) q[3];
sx q[3];
rz(-2.4628277) q[3];
sx q[3];
rz(0.15798727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79992646) q[2];
sx q[2];
rz(-3.0604) q[2];
sx q[2];
rz(-2.9426835) q[2];
rz(0.79814664) q[3];
sx q[3];
rz(-1.9559559) q[3];
sx q[3];
rz(-1.859349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9195008) q[0];
sx q[0];
rz(-1.5380479) q[0];
sx q[0];
rz(2.0684239) q[0];
rz(2.3566698) q[1];
sx q[1];
rz(-2.2833318) q[1];
sx q[1];
rz(-2.2699184) q[1];
rz(1.9368108) q[2];
sx q[2];
rz(-3.0095965) q[2];
sx q[2];
rz(-2.1684564) q[2];
rz(-0.34134117) q[3];
sx q[3];
rz(-0.69820709) q[3];
sx q[3];
rz(0.38521614) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
