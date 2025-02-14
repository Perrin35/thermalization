OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1342993) q[0];
sx q[0];
rz(-1.1987885) q[0];
sx q[0];
rz(2.3564136) q[0];
rz(2.5685318) q[1];
sx q[1];
rz(-0.95102024) q[1];
sx q[1];
rz(-2.5101488) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9438659) q[0];
sx q[0];
rz(-2.336204) q[0];
sx q[0];
rz(0.95279537) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72452568) q[2];
sx q[2];
rz(-0.93225485) q[2];
sx q[2];
rz(2.2341626) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3062485) q[1];
sx q[1];
rz(-2.1901096) q[1];
sx q[1];
rz(2.1080416) q[1];
x q[2];
rz(-1.6699766) q[3];
sx q[3];
rz(-2.0137312) q[3];
sx q[3];
rz(0.25275005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.50420061) q[2];
sx q[2];
rz(-1.4376419) q[2];
sx q[2];
rz(0.27153095) q[2];
rz(3.0951989) q[3];
sx q[3];
rz(-1.4103187) q[3];
sx q[3];
rz(2.7739649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72535998) q[0];
sx q[0];
rz(-3.1144436) q[0];
sx q[0];
rz(2.695988) q[0];
rz(-2.7815869) q[1];
sx q[1];
rz(-1.5574734) q[1];
sx q[1];
rz(0.97703591) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1030658) q[0];
sx q[0];
rz(-1.0630106) q[0];
sx q[0];
rz(-0.610791) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4233098) q[2];
sx q[2];
rz(-0.86330763) q[2];
sx q[2];
rz(-1.8607418) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9865446) q[1];
sx q[1];
rz(-0.98216559) q[1];
sx q[1];
rz(-0.18208336) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3958443) q[3];
sx q[3];
rz(-0.85674113) q[3];
sx q[3];
rz(-0.24859771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6147324) q[2];
sx q[2];
rz(-1.7975668) q[2];
sx q[2];
rz(-0.942222) q[2];
rz(2.6451536) q[3];
sx q[3];
rz(-1.6783291) q[3];
sx q[3];
rz(1.0676603) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0838858) q[0];
sx q[0];
rz(-1.3297465) q[0];
sx q[0];
rz(0.65621334) q[0];
rz(0.48037848) q[1];
sx q[1];
rz(-2.1800133) q[1];
sx q[1];
rz(0.61294714) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543535) q[0];
sx q[0];
rz(-1.6677852) q[0];
sx q[0];
rz(1.9796275) q[0];
rz(-0.0080994227) q[2];
sx q[2];
rz(-0.73609422) q[2];
sx q[2];
rz(-1.8437633) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7483116) q[1];
sx q[1];
rz(-1.5419811) q[1];
sx q[1];
rz(-0.033571413) q[1];
rz(1.5016236) q[3];
sx q[3];
rz(-1.072848) q[3];
sx q[3];
rz(2.43735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2059242) q[2];
sx q[2];
rz(-1.0574477) q[2];
sx q[2];
rz(0.52440468) q[2];
rz(-2.2762401) q[3];
sx q[3];
rz(-1.0582346) q[3];
sx q[3];
rz(-2.481148) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.653729) q[0];
sx q[0];
rz(-1.3069343) q[0];
sx q[0];
rz(-3.0453239) q[0];
rz(3.1277711) q[1];
sx q[1];
rz(-2.6410069) q[1];
sx q[1];
rz(1.0523419) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5139363) q[0];
sx q[0];
rz(-0.13051438) q[0];
sx q[0];
rz(1.9399727) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8554535) q[2];
sx q[2];
rz(-2.1026662) q[2];
sx q[2];
rz(-1.4653613) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2675002) q[1];
sx q[1];
rz(-0.37768957) q[1];
sx q[1];
rz(0.51376359) q[1];
rz(-pi) q[2];
rz(-0.035798612) q[3];
sx q[3];
rz(-1.0455264) q[3];
sx q[3];
rz(-2.4095499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.72157613) q[2];
sx q[2];
rz(-2.6105011) q[2];
sx q[2];
rz(0.51900035) q[2];
rz(1.0985993) q[3];
sx q[3];
rz(-1.6853305) q[3];
sx q[3];
rz(0.73002735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62243432) q[0];
sx q[0];
rz(-0.19304481) q[0];
sx q[0];
rz(-2.4073507) q[0];
rz(1.2892067) q[1];
sx q[1];
rz(-2.0964918) q[1];
sx q[1];
rz(2.2330914) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8216128) q[0];
sx q[0];
rz(-2.2083518) q[0];
sx q[0];
rz(2.430401) q[0];
rz(-pi) q[1];
rz(1.7041358) q[2];
sx q[2];
rz(-2.1380084) q[2];
sx q[2];
rz(0.043616488) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.31181609) q[1];
sx q[1];
rz(-0.83925345) q[1];
sx q[1];
rz(-3.0481382) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8374422) q[3];
sx q[3];
rz(-1.5943267) q[3];
sx q[3];
rz(-1.2465828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.044540731) q[2];
sx q[2];
rz(-2.3553607) q[2];
sx q[2];
rz(2.457974) q[2];
rz(2.0475856) q[3];
sx q[3];
rz(-1.3235599) q[3];
sx q[3];
rz(1.5723642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053726824) q[0];
sx q[0];
rz(-1.6931067) q[0];
sx q[0];
rz(-2.7159765) q[0];
rz(2.5348237) q[1];
sx q[1];
rz(-0.68777045) q[1];
sx q[1];
rz(-1.2026131) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86798651) q[0];
sx q[0];
rz(-0.79279723) q[0];
sx q[0];
rz(-2.2895564) q[0];
x q[1];
rz(-3.075766) q[2];
sx q[2];
rz(-1.5562686) q[2];
sx q[2];
rz(-0.17958497) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7640671) q[1];
sx q[1];
rz(-2.0646618) q[1];
sx q[1];
rz(-0.74158828) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7814066) q[3];
sx q[3];
rz(-2.8202116) q[3];
sx q[3];
rz(1.7226294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6732424) q[2];
sx q[2];
rz(-0.52933669) q[2];
sx q[2];
rz(2.5128515) q[2];
rz(-0.33744234) q[3];
sx q[3];
rz(-1.7818261) q[3];
sx q[3];
rz(2.3745645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5479946) q[0];
sx q[0];
rz(-3.0715521) q[0];
sx q[0];
rz(1.8016169) q[0];
rz(-3.0361259) q[1];
sx q[1];
rz(-2.1601845) q[1];
sx q[1];
rz(2.9973105) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54050416) q[0];
sx q[0];
rz(-1.2799037) q[0];
sx q[0];
rz(2.7379047) q[0];
rz(-pi) q[1];
rz(2.9469523) q[2];
sx q[2];
rz(-0.87207672) q[2];
sx q[2];
rz(-0.34572225) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6444417) q[1];
sx q[1];
rz(-1.9195917) q[1];
sx q[1];
rz(1.3494874) q[1];
x q[2];
rz(2.1142247) q[3];
sx q[3];
rz(-0.97350922) q[3];
sx q[3];
rz(-1.3940108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0353126) q[2];
sx q[2];
rz(-1.148618) q[2];
sx q[2];
rz(2.3036352) q[2];
rz(0.49977866) q[3];
sx q[3];
rz(-0.79334799) q[3];
sx q[3];
rz(-2.4404073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4226294) q[0];
sx q[0];
rz(-2.0327649) q[0];
sx q[0];
rz(-0.75123373) q[0];
rz(-1.2535837) q[1];
sx q[1];
rz(-1.8389849) q[1];
sx q[1];
rz(1.368103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73117646) q[0];
sx q[0];
rz(-2.3623657) q[0];
sx q[0];
rz(2.221629) q[0];
x q[1];
rz(1.8644731) q[2];
sx q[2];
rz(-1.8741359) q[2];
sx q[2];
rz(1.9544698) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.320619) q[1];
sx q[1];
rz(-1.486196) q[1];
sx q[1];
rz(2.3449347) q[1];
rz(-pi) q[2];
rz(-2.1001216) q[3];
sx q[3];
rz(-2.4987037) q[3];
sx q[3];
rz(2.3704308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3086243) q[2];
sx q[2];
rz(-0.37165752) q[2];
sx q[2];
rz(-0.0011681636) q[2];
rz(0.33231169) q[3];
sx q[3];
rz(-1.6806335) q[3];
sx q[3];
rz(-2.6619679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7697656) q[0];
sx q[0];
rz(-1.9244939) q[0];
sx q[0];
rz(-1.0960854) q[0];
rz(-1.7425884) q[1];
sx q[1];
rz(-1.5053446) q[1];
sx q[1];
rz(2.2217264) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4324661) q[0];
sx q[0];
rz(-1.5282915) q[0];
sx q[0];
rz(2.0086238) q[0];
rz(-pi) q[1];
rz(0.28025744) q[2];
sx q[2];
rz(-1.9748903) q[2];
sx q[2];
rz(2.0660482) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.49860543) q[1];
sx q[1];
rz(-2.4772493) q[1];
sx q[1];
rz(-0.24533116) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3289973) q[3];
sx q[3];
rz(-2.5957172) q[3];
sx q[3];
rz(0.98009566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2190735) q[2];
sx q[2];
rz(-2.8963431) q[2];
sx q[2];
rz(1.5101439) q[2];
rz(-0.63589969) q[3];
sx q[3];
rz(-0.70521516) q[3];
sx q[3];
rz(-0.51457921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9566327) q[0];
sx q[0];
rz(-2.2745467) q[0];
sx q[0];
rz(1.4956723) q[0];
rz(3.0236859) q[1];
sx q[1];
rz(-1.7472183) q[1];
sx q[1];
rz(-2.8848021) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2599597) q[0];
sx q[0];
rz(-1.6266394) q[0];
sx q[0];
rz(2.0042438) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58509201) q[2];
sx q[2];
rz(-1.3199542) q[2];
sx q[2];
rz(-2.5426007) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6424055) q[1];
sx q[1];
rz(-0.39246628) q[1];
sx q[1];
rz(-2.3981903) q[1];
rz(-pi) q[2];
x q[2];
rz(2.58415) q[3];
sx q[3];
rz(-2.5187105) q[3];
sx q[3];
rz(1.0475782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7439277) q[2];
sx q[2];
rz(-0.74349082) q[2];
sx q[2];
rz(-0.94997722) q[2];
rz(2.5395565) q[3];
sx q[3];
rz(-1.2369912) q[3];
sx q[3];
rz(-0.33815798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6947185) q[0];
sx q[0];
rz(-2.5304486) q[0];
sx q[0];
rz(1.9421938) q[0];
rz(-2.1318204) q[1];
sx q[1];
rz(-0.92242454) q[1];
sx q[1];
rz(-0.19933272) q[1];
rz(-1.2620325) q[2];
sx q[2];
rz(-1.8539351) q[2];
sx q[2];
rz(-1.2507982) q[2];
rz(0.39691858) q[3];
sx q[3];
rz(-0.81868681) q[3];
sx q[3];
rz(0.80358748) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
