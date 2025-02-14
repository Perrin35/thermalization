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
rz(1.4242564) q[0];
sx q[0];
rz(-1.2196701) q[0];
sx q[0];
rz(-0.57504672) q[0];
rz(-2.9991034) q[1];
sx q[1];
rz(-1.2187076) q[1];
sx q[1];
rz(0.86427468) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48102114) q[0];
sx q[0];
rz(-0.29566524) q[0];
sx q[0];
rz(-0.85938248) q[0];
rz(0.15344413) q[2];
sx q[2];
rz(-2.1965003) q[2];
sx q[2];
rz(0.99391711) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4839403) q[1];
sx q[1];
rz(-0.53002702) q[1];
sx q[1];
rz(0.0040667314) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0686758) q[3];
sx q[3];
rz(-2.0808633) q[3];
sx q[3];
rz(-0.74397457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.064726) q[2];
sx q[2];
rz(-1.4996424) q[2];
sx q[2];
rz(1.0092674) q[2];
rz(0.81389728) q[3];
sx q[3];
rz(-2.8129357) q[3];
sx q[3];
rz(2.8024659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086804248) q[0];
sx q[0];
rz(-1.7330994) q[0];
sx q[0];
rz(2.8993697) q[0];
rz(1.9107266) q[1];
sx q[1];
rz(-0.63175285) q[1];
sx q[1];
rz(2.3435074) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.674192) q[0];
sx q[0];
rz(-2.5134006) q[0];
sx q[0];
rz(-0.19285658) q[0];
rz(-pi) q[1];
rz(1.2276137) q[2];
sx q[2];
rz(-1.7084885) q[2];
sx q[2];
rz(-2.340449) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6769413) q[1];
sx q[1];
rz(-1.2124227) q[1];
sx q[1];
rz(-0.48848571) q[1];
x q[2];
rz(0.40753461) q[3];
sx q[3];
rz(-0.38713249) q[3];
sx q[3];
rz(-2.3705496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65608394) q[2];
sx q[2];
rz(-1.069243) q[2];
sx q[2];
rz(-2.3089224) q[2];
rz(-0.63140702) q[3];
sx q[3];
rz(-0.74149817) q[3];
sx q[3];
rz(3.1413063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0176508) q[0];
sx q[0];
rz(-2.232382) q[0];
sx q[0];
rz(0.18774524) q[0];
rz(-2.2979157) q[1];
sx q[1];
rz(-1.2742821) q[1];
sx q[1];
rz(0.63293308) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3819861) q[0];
sx q[0];
rz(-1.6416024) q[0];
sx q[0];
rz(-1.9340865) q[0];
rz(1.2212513) q[2];
sx q[2];
rz(-2.2270906) q[2];
sx q[2];
rz(1.2906769) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.17281118) q[1];
sx q[1];
rz(-2.079614) q[1];
sx q[1];
rz(-1.225586) q[1];
rz(-pi) q[2];
rz(0.12944451) q[3];
sx q[3];
rz(-1.9151546) q[3];
sx q[3];
rz(0.5705006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26744276) q[2];
sx q[2];
rz(-0.94620693) q[2];
sx q[2];
rz(1.2713185) q[2];
rz(1.4960131) q[3];
sx q[3];
rz(-1.6451903) q[3];
sx q[3];
rz(-0.97150826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2520168) q[0];
sx q[0];
rz(-2.3887964) q[0];
sx q[0];
rz(0.65943199) q[0];
rz(1.8239498) q[1];
sx q[1];
rz(-1.3392071) q[1];
sx q[1];
rz(2.3514294) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5555429) q[0];
sx q[0];
rz(-1.2359185) q[0];
sx q[0];
rz(1.4396776) q[0];
x q[1];
rz(-2.1151562) q[2];
sx q[2];
rz(-3.1033278) q[2];
sx q[2];
rz(-2.1452034) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9875264) q[1];
sx q[1];
rz(-0.48168698) q[1];
sx q[1];
rz(0.14194004) q[1];
x q[2];
rz(0.55410093) q[3];
sx q[3];
rz(-1.3375207) q[3];
sx q[3];
rz(-1.8221888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.059375199) q[2];
sx q[2];
rz(-1.6148115) q[2];
sx q[2];
rz(0.29602948) q[2];
rz(2.5791903) q[3];
sx q[3];
rz(-0.86807576) q[3];
sx q[3];
rz(-3.0276022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89428467) q[0];
sx q[0];
rz(-1.9534651) q[0];
sx q[0];
rz(-2.9344015) q[0];
rz(1.0317135) q[1];
sx q[1];
rz(-1.9887911) q[1];
sx q[1];
rz(-1.656104) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0061917) q[0];
sx q[0];
rz(-1.1814607) q[0];
sx q[0];
rz(-0.72442309) q[0];
rz(-pi) q[1];
rz(-2.1319415) q[2];
sx q[2];
rz(-1.3232987) q[2];
sx q[2];
rz(-1.1010608) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7149026) q[1];
sx q[1];
rz(-1.443092) q[1];
sx q[1];
rz(2.7903778) q[1];
rz(0.84019227) q[3];
sx q[3];
rz(-1.4903791) q[3];
sx q[3];
rz(-1.059746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1032054) q[2];
sx q[2];
rz(-0.81192553) q[2];
sx q[2];
rz(-1.1078328) q[2];
rz(2.2191018) q[3];
sx q[3];
rz(-1.731571) q[3];
sx q[3];
rz(-2.8973575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.51782) q[0];
sx q[0];
rz(-0.51384846) q[0];
sx q[0];
rz(-2.9129831) q[0];
rz(-0.27433431) q[1];
sx q[1];
rz(-2.3286596) q[1];
sx q[1];
rz(2.6913604) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58670973) q[0];
sx q[0];
rz(-2.3392117) q[0];
sx q[0];
rz(-0.35361617) q[0];
x q[1];
rz(-1.8000431) q[2];
sx q[2];
rz(-1.5783797) q[2];
sx q[2];
rz(0.46310616) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6630312) q[1];
sx q[1];
rz(-1.3197462) q[1];
sx q[1];
rz(0.51086225) q[1];
rz(-pi) q[2];
rz(2.130533) q[3];
sx q[3];
rz(-0.478906) q[3];
sx q[3];
rz(-0.68747179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7776514) q[2];
sx q[2];
rz(-0.55299091) q[2];
sx q[2];
rz(3.1138163) q[2];
rz(2.3751496) q[3];
sx q[3];
rz(-1.1602297) q[3];
sx q[3];
rz(0.69507039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.22208333) q[0];
sx q[0];
rz(-1.6351901) q[0];
sx q[0];
rz(-0.62244225) q[0];
rz(2.4355603) q[1];
sx q[1];
rz(-2.249554) q[1];
sx q[1];
rz(0.58696729) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7919803) q[0];
sx q[0];
rz(-1.1328348) q[0];
sx q[0];
rz(2.1823606) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48611792) q[2];
sx q[2];
rz(-1.3335506) q[2];
sx q[2];
rz(-3.089774) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.5712292) q[1];
sx q[1];
rz(-1.0652568) q[1];
sx q[1];
rz(-1.4897219) q[1];
rz(-1.582828) q[3];
sx q[3];
rz(-2.4295252) q[3];
sx q[3];
rz(-2.4648417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.065757699) q[2];
sx q[2];
rz(-1.4328052) q[2];
sx q[2];
rz(2.5471121) q[2];
rz(-2.8822656) q[3];
sx q[3];
rz(-1.2649053) q[3];
sx q[3];
rz(1.2172786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1031621) q[0];
sx q[0];
rz(-2.2362464) q[0];
sx q[0];
rz(-2.5627947) q[0];
rz(-1.831306) q[1];
sx q[1];
rz(-1.2625887) q[1];
sx q[1];
rz(-2.328918) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2437163) q[0];
sx q[0];
rz(-1.0930632) q[0];
sx q[0];
rz(-1.7327139) q[0];
rz(1.0092998) q[2];
sx q[2];
rz(-1.7914869) q[2];
sx q[2];
rz(-2.7796628) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1393237) q[1];
sx q[1];
rz(-1.9319169) q[1];
sx q[1];
rz(0.72092542) q[1];
rz(0.854354) q[3];
sx q[3];
rz(-2.4378447) q[3];
sx q[3];
rz(3.0643413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3506713) q[2];
sx q[2];
rz(-1.8427883) q[2];
sx q[2];
rz(2.217963) q[2];
rz(-2.1212497) q[3];
sx q[3];
rz(-0.89710051) q[3];
sx q[3];
rz(-3.0237107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.0608805) q[0];
sx q[0];
rz(-1.7055644) q[0];
sx q[0];
rz(-1.4870148) q[0];
rz(3.0912073) q[1];
sx q[1];
rz(-0.81704187) q[1];
sx q[1];
rz(2.494716) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4934702) q[0];
sx q[0];
rz(-0.97705829) q[0];
sx q[0];
rz(0.52072449) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5570103) q[2];
sx q[2];
rz(-2.3858831) q[2];
sx q[2];
rz(-0.85697714) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8343492) q[1];
sx q[1];
rz(-1.3615047) q[1];
sx q[1];
rz(0.41469519) q[1];
x q[2];
rz(1.3573801) q[3];
sx q[3];
rz(-1.4239256) q[3];
sx q[3];
rz(-1.980912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6369624) q[2];
sx q[2];
rz(-1.5550104) q[2];
sx q[2];
rz(-0.75616765) q[2];
rz(-1.470083) q[3];
sx q[3];
rz(-2.3556637) q[3];
sx q[3];
rz(0.80529958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1038372) q[0];
sx q[0];
rz(-1.8917731) q[0];
sx q[0];
rz(-2.0334429) q[0];
rz(-2.8773819) q[1];
sx q[1];
rz(-1.4690396) q[1];
sx q[1];
rz(0.60233751) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8648351) q[0];
sx q[0];
rz(-1.2981155) q[0];
sx q[0];
rz(-2.9972267) q[0];
rz(-pi) q[1];
rz(-0.71508821) q[2];
sx q[2];
rz(-0.23221604) q[2];
sx q[2];
rz(0.84399904) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70302898) q[1];
sx q[1];
rz(-1.6191585) q[1];
sx q[1];
rz(-0.88367422) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65091316) q[3];
sx q[3];
rz(-1.7554566) q[3];
sx q[3];
rz(2.1391275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.63310528) q[2];
sx q[2];
rz(-2.5999531) q[2];
sx q[2];
rz(2.2701021) q[2];
rz(1.9320711) q[3];
sx q[3];
rz(-0.28035527) q[3];
sx q[3];
rz(2.5476294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4393944) q[0];
sx q[0];
rz(-1.7560503) q[0];
sx q[0];
rz(-0.51881292) q[0];
rz(2.5595472) q[1];
sx q[1];
rz(-1.1036292) q[1];
sx q[1];
rz(-1.6694952) q[1];
rz(-1.9459361) q[2];
sx q[2];
rz(-0.87004253) q[2];
sx q[2];
rz(-0.042849356) q[2];
rz(2.9357435) q[3];
sx q[3];
rz(-2.4424853) q[3];
sx q[3];
rz(1.427099) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
