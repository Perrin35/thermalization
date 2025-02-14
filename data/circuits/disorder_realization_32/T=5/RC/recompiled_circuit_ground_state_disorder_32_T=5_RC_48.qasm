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
rz(-0.78517908) q[0];
rz(-0.5730609) q[1];
sx q[1];
rz(-2.1905724) q[1];
sx q[1];
rz(-0.63144389) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60053958) q[0];
sx q[0];
rz(-2.1990393) q[0];
sx q[0];
rz(2.598935) q[0];
rz(-pi) q[1];
rz(-0.84191834) q[2];
sx q[2];
rz(-0.92570423) q[2];
sx q[2];
rz(-3.0708714) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.50791477) q[1];
sx q[1];
rz(-0.79601015) q[1];
sx q[1];
rz(-2.518954) q[1];
x q[2];
rz(-0.20578014) q[3];
sx q[3];
rz(-0.45318402) q[3];
sx q[3];
rz(-3.1169716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.637392) q[2];
sx q[2];
rz(-1.4376419) q[2];
sx q[2];
rz(-2.8700617) q[2];
rz(3.0951989) q[3];
sx q[3];
rz(-1.4103187) q[3];
sx q[3];
rz(-0.36762777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72535998) q[0];
sx q[0];
rz(-0.027149057) q[0];
sx q[0];
rz(-0.44560462) q[0];
rz(-0.36000571) q[1];
sx q[1];
rz(-1.5841192) q[1];
sx q[1];
rz(-2.1645567) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038526857) q[0];
sx q[0];
rz(-2.0785821) q[0];
sx q[0];
rz(-2.5308017) q[0];
x q[1];
rz(2.4196834) q[2];
sx q[2];
rz(-1.0470265) q[2];
sx q[2];
rz(0.80654752) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9865446) q[1];
sx q[1];
rz(-0.98216559) q[1];
sx q[1];
rz(-2.9595093) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70313022) q[3];
sx q[3];
rz(-2.1092013) q[3];
sx q[3];
rz(2.3634274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.52686024) q[2];
sx q[2];
rz(-1.7975668) q[2];
sx q[2];
rz(-0.942222) q[2];
rz(0.49643907) q[3];
sx q[3];
rz(-1.4632635) q[3];
sx q[3];
rz(1.0676603) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0577069) q[0];
sx q[0];
rz(-1.8118462) q[0];
sx q[0];
rz(-0.65621334) q[0];
rz(2.6612142) q[1];
sx q[1];
rz(-0.96157938) q[1];
sx q[1];
rz(0.61294714) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1663294) q[0];
sx q[0];
rz(-2.722046) q[0];
sx q[0];
rz(1.3307721) q[0];
x q[1];
rz(-3.1334932) q[2];
sx q[2];
rz(-2.4054984) q[2];
sx q[2];
rz(1.2978293) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9631098) q[1];
sx q[1];
rz(-1.5372389) q[1];
sx q[1];
rz(-1.5996278) q[1];
x q[2];
rz(-1.6399691) q[3];
sx q[3];
rz(-1.072848) q[3];
sx q[3];
rz(-0.70424265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9356685) q[2];
sx q[2];
rz(-2.084145) q[2];
sx q[2];
rz(-2.617188) q[2];
rz(-2.2762401) q[3];
sx q[3];
rz(-1.0582346) q[3];
sx q[3];
rz(-2.481148) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48786369) q[0];
sx q[0];
rz(-1.3069343) q[0];
sx q[0];
rz(3.0453239) q[0];
rz(-3.1277711) q[1];
sx q[1];
rz(-0.50058573) q[1];
sx q[1];
rz(-2.0892508) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5139363) q[0];
sx q[0];
rz(-0.13051438) q[0];
sx q[0];
rz(-1.9399727) q[0];
x q[1];
rz(-1.123548) q[2];
sx q[2];
rz(-2.5442225) q[2];
sx q[2];
rz(-1.9910461) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7218771) q[1];
sx q[1];
rz(-1.8977563) q[1];
sx q[1];
rz(1.7633597) q[1];
rz(-pi) q[2];
rz(-3.105794) q[3];
sx q[3];
rz(-1.0455264) q[3];
sx q[3];
rz(-0.73204277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4200165) q[2];
sx q[2];
rz(-0.53109157) q[2];
sx q[2];
rz(-2.6225923) q[2];
rz(-1.0985993) q[3];
sx q[3];
rz(-1.6853305) q[3];
sx q[3];
rz(-0.73002735) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5191583) q[0];
sx q[0];
rz(-0.19304481) q[0];
sx q[0];
rz(0.73424196) q[0];
rz(1.8523859) q[1];
sx q[1];
rz(-1.0451008) q[1];
sx q[1];
rz(-0.9085013) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4168981) q[0];
sx q[0];
rz(-1.0186581) q[0];
sx q[0];
rz(0.79663251) q[0];
x q[1];
rz(1.4374569) q[2];
sx q[2];
rz(-2.1380084) q[2];
sx q[2];
rz(3.0979762) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.1724194) q[1];
sx q[1];
rz(-2.4052019) q[1];
sx q[1];
rz(1.6743771) q[1];
x q[2];
rz(-3.1172006) q[3];
sx q[3];
rz(-1.304226) q[3];
sx q[3];
rz(-0.3177869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0970519) q[2];
sx q[2];
rz(-2.3553607) q[2];
sx q[2];
rz(-0.68361863) q[2];
rz(1.094007) q[3];
sx q[3];
rz(-1.3235599) q[3];
sx q[3];
rz(1.5692284) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053726824) q[0];
sx q[0];
rz(-1.448486) q[0];
sx q[0];
rz(-0.42561612) q[0];
rz(-0.60676891) q[1];
sx q[1];
rz(-0.68777045) q[1];
sx q[1];
rz(1.9389796) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9894597) q[0];
sx q[0];
rz(-1.0826064) q[0];
sx q[0];
rz(2.2230985) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5853556) q[2];
sx q[2];
rz(-1.636616) q[2];
sx q[2];
rz(1.3902537) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9252117) q[1];
sx q[1];
rz(-2.2077476) q[1];
sx q[1];
rz(-0.94016192) q[1];
rz(-pi) q[2];
rz(2.909502) q[3];
sx q[3];
rz(-1.795139) q[3];
sx q[3];
rz(2.2346131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6732424) q[2];
sx q[2];
rz(-2.612256) q[2];
sx q[2];
rz(-2.5128515) q[2];
rz(0.33744234) q[3];
sx q[3];
rz(-1.3597666) q[3];
sx q[3];
rz(-0.76702816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5479946) q[0];
sx q[0];
rz(-0.070040528) q[0];
sx q[0];
rz(-1.8016169) q[0];
rz(-3.0361259) q[1];
sx q[1];
rz(-0.98140812) q[1];
sx q[1];
rz(-2.9973105) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54050416) q[0];
sx q[0];
rz(-1.2799037) q[0];
sx q[0];
rz(2.7379047) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86267306) q[2];
sx q[2];
rz(-1.7194334) q[2];
sx q[2];
rz(1.3512063) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.99690106) q[1];
sx q[1];
rz(-1.7785774) q[1];
sx q[1];
rz(-2.7847915) q[1];
x q[2];
rz(2.470131) q[3];
sx q[3];
rz(-2.0125768) q[3];
sx q[3];
rz(-2.9908671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0353126) q[2];
sx q[2];
rz(-1.148618) q[2];
sx q[2];
rz(-2.3036352) q[2];
rz(2.641814) q[3];
sx q[3];
rz(-2.3482447) q[3];
sx q[3];
rz(0.70118538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4226294) q[0];
sx q[0];
rz(-2.0327649) q[0];
sx q[0];
rz(2.3903589) q[0];
rz(1.8880089) q[1];
sx q[1];
rz(-1.3026078) q[1];
sx q[1];
rz(-1.368103) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7984895) q[0];
sx q[0];
rz(-2.0105848) q[0];
sx q[0];
rz(-2.2368311) q[0];
rz(0.31603916) q[2];
sx q[2];
rz(-1.2908986) q[2];
sx q[2];
rz(-2.6678277) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66361278) q[1];
sx q[1];
rz(-2.3638032) q[1];
sx q[1];
rz(1.4500834) q[1];
x q[2];
rz(-2.1447324) q[3];
sx q[3];
rz(-1.2632476) q[3];
sx q[3];
rz(2.779863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.83296835) q[2];
sx q[2];
rz(-0.37165752) q[2];
sx q[2];
rz(3.1404245) q[2];
rz(-0.33231169) q[3];
sx q[3];
rz(-1.4609591) q[3];
sx q[3];
rz(0.47962475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.371827) q[0];
sx q[0];
rz(-1.2170987) q[0];
sx q[0];
rz(-2.0455072) q[0];
rz(1.3990043) q[1];
sx q[1];
rz(-1.636248) q[1];
sx q[1];
rz(-2.2217264) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8815589) q[0];
sx q[0];
rz(-1.1333916) q[0];
sx q[0];
rz(-0.046925515) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1521173) q[2];
sx q[2];
rz(-1.3136465) q[2];
sx q[2];
rz(2.5336483) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.49860543) q[1];
sx q[1];
rz(-2.4772493) q[1];
sx q[1];
rz(-2.8962615) q[1];
x q[2];
rz(2.1036622) q[3];
sx q[3];
rz(-1.4461596) q[3];
sx q[3];
rz(-0.79844284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9225191) q[2];
sx q[2];
rz(-2.8963431) q[2];
sx q[2];
rz(-1.5101439) q[2];
rz(-0.63589969) q[3];
sx q[3];
rz(-2.4363775) q[3];
sx q[3];
rz(0.51457921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9566327) q[0];
sx q[0];
rz(-2.2745467) q[0];
sx q[0];
rz(1.6459203) q[0];
rz(3.0236859) q[1];
sx q[1];
rz(-1.7472183) q[1];
sx q[1];
rz(-2.8848021) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3366616) q[0];
sx q[0];
rz(-1.1380702) q[0];
sx q[0];
rz(-0.061519775) q[0];
rz(2.5565006) q[2];
sx q[2];
rz(-1.8216385) q[2];
sx q[2];
rz(-0.59899194) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8578882) q[1];
sx q[1];
rz(-1.285375) q[1];
sx q[1];
rz(1.2976451) q[1];
x q[2];
rz(-1.9339233) q[3];
sx q[3];
rz(-2.0887017) q[3];
sx q[3];
rz(-0.3929485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7439277) q[2];
sx q[2];
rz(-2.3981018) q[2];
sx q[2];
rz(0.94997722) q[2];
rz(0.60203612) q[3];
sx q[3];
rz(-1.9046015) q[3];
sx q[3];
rz(2.8034347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6947185) q[0];
sx q[0];
rz(-0.61114408) q[0];
sx q[0];
rz(-1.1993988) q[0];
rz(-2.1318204) q[1];
sx q[1];
rz(-0.92242454) q[1];
sx q[1];
rz(-0.19933272) q[1];
rz(0.29640167) q[2];
sx q[2];
rz(-1.2747073) q[2];
sx q[2];
rz(-2.7327197) q[2];
rz(-1.9626403) q[3];
sx q[3];
rz(-0.83189356) q[3];
sx q[3];
rz(-2.8883285) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
