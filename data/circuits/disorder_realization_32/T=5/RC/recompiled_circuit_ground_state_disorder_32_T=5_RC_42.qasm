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
rz(2.1905724) q[1];
sx q[1];
rz(11.934927) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60053958) q[0];
sx q[0];
rz(-0.94255336) q[0];
sx q[0];
rz(0.54265768) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2996743) q[2];
sx q[2];
rz(-2.2158884) q[2];
sx q[2];
rz(-3.0708714) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8353442) q[1];
sx q[1];
rz(-2.1901096) q[1];
sx q[1];
rz(-2.1080416) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6699766) q[3];
sx q[3];
rz(-2.0137312) q[3];
sx q[3];
rz(2.8888426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.637392) q[2];
sx q[2];
rz(-1.7039508) q[2];
sx q[2];
rz(-2.8700617) q[2];
rz(-3.0951989) q[3];
sx q[3];
rz(-1.4103187) q[3];
sx q[3];
rz(0.36762777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4162327) q[0];
sx q[0];
rz(-3.1144436) q[0];
sx q[0];
rz(2.695988) q[0];
rz(-2.7815869) q[1];
sx q[1];
rz(-1.5574734) q[1];
sx q[1];
rz(0.97703591) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8603823) q[0];
sx q[0];
rz(-2.095725) q[0];
sx q[0];
rz(0.97411823) q[0];
x q[1];
rz(-0.71828281) q[2];
sx q[2];
rz(-2.278285) q[2];
sx q[2];
rz(1.2808509) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9865446) q[1];
sx q[1];
rz(-2.1594271) q[1];
sx q[1];
rz(0.18208336) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74574836) q[3];
sx q[3];
rz(-2.2848515) q[3];
sx q[3];
rz(0.24859771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6147324) q[2];
sx q[2];
rz(-1.7975668) q[2];
sx q[2];
rz(-2.1993707) q[2];
rz(-2.6451536) q[3];
sx q[3];
rz(-1.6783291) q[3];
sx q[3];
rz(-1.0676603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0838858) q[0];
sx q[0];
rz(-1.3297465) q[0];
sx q[0];
rz(-0.65621334) q[0];
rz(2.6612142) q[1];
sx q[1];
rz(-2.1800133) q[1];
sx q[1];
rz(2.5286455) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71343173) q[0];
sx q[0];
rz(-1.9775922) q[0];
sx q[0];
rz(3.0359546) q[0];
x q[1];
rz(1.5634588) q[2];
sx q[2];
rz(-0.83473182) q[2];
sx q[2];
rz(1.8328345) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7483116) q[1];
sx q[1];
rz(-1.5419811) q[1];
sx q[1];
rz(0.033571413) q[1];
x q[2];
rz(-0.49895371) q[3];
sx q[3];
rz(-1.5100347) q[3];
sx q[3];
rz(2.2419597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9356685) q[2];
sx q[2];
rz(-2.084145) q[2];
sx q[2];
rz(-2.617188) q[2];
rz(2.2762401) q[3];
sx q[3];
rz(-2.083358) q[3];
sx q[3];
rz(0.66044468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48786369) q[0];
sx q[0];
rz(-1.8346584) q[0];
sx q[0];
rz(-0.096268795) q[0];
rz(-0.013821566) q[1];
sx q[1];
rz(-0.50058573) q[1];
sx q[1];
rz(2.0892508) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6276564) q[0];
sx q[0];
rz(-0.13051438) q[0];
sx q[0];
rz(-1.2016199) q[0];
rz(2.8554535) q[2];
sx q[2];
rz(-1.0389265) q[2];
sx q[2];
rz(-1.6762313) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7218771) q[1];
sx q[1];
rz(-1.8977563) q[1];
sx q[1];
rz(1.378233) q[1];
rz(-pi) q[2];
rz(-1.0452483) q[3];
sx q[3];
rz(-1.6017672) q[3];
sx q[3];
rz(-2.3207959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4200165) q[2];
sx q[2];
rz(-0.53109157) q[2];
sx q[2];
rz(-0.51900035) q[2];
rz(-1.0985993) q[3];
sx q[3];
rz(-1.6853305) q[3];
sx q[3];
rz(-0.73002735) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62243432) q[0];
sx q[0];
rz(-2.9485478) q[0];
sx q[0];
rz(-2.4073507) q[0];
rz(-1.2892067) q[1];
sx q[1];
rz(-2.0964918) q[1];
sx q[1];
rz(-2.2330914) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72469456) q[0];
sx q[0];
rz(-1.0186581) q[0];
sx q[0];
rz(-0.79663251) q[0];
x q[1];
rz(-1.7041358) q[2];
sx q[2];
rz(-1.0035842) q[2];
sx q[2];
rz(0.043616488) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.31181609) q[1];
sx q[1];
rz(-2.3023392) q[1];
sx q[1];
rz(-0.093454425) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1172006) q[3];
sx q[3];
rz(-1.304226) q[3];
sx q[3];
rz(-2.8238058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.044540731) q[2];
sx q[2];
rz(-2.3553607) q[2];
sx q[2];
rz(-0.68361863) q[2];
rz(-1.094007) q[3];
sx q[3];
rz(-1.3235599) q[3];
sx q[3];
rz(-1.5692284) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0878658) q[0];
sx q[0];
rz(-1.6931067) q[0];
sx q[0];
rz(-0.42561612) q[0];
rz(2.5348237) q[1];
sx q[1];
rz(-2.4538222) q[1];
sx q[1];
rz(1.2026131) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9894597) q[0];
sx q[0];
rz(-1.0826064) q[0];
sx q[0];
rz(0.91849416) q[0];
x q[1];
rz(-1.5853556) q[2];
sx q[2];
rz(-1.5049767) q[2];
sx q[2];
rz(1.3902537) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.21638097) q[1];
sx q[1];
rz(-0.93384508) q[1];
sx q[1];
rz(2.2014307) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.909502) q[3];
sx q[3];
rz(-1.3464536) q[3];
sx q[3];
rz(-0.90697955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6732424) q[2];
sx q[2];
rz(-0.52933669) q[2];
sx q[2];
rz(-2.5128515) q[2];
rz(-2.8041503) q[3];
sx q[3];
rz(-1.7818261) q[3];
sx q[3];
rz(-2.3745645) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5479946) q[0];
sx q[0];
rz(-0.070040528) q[0];
sx q[0];
rz(1.8016169) q[0];
rz(-0.10546671) q[1];
sx q[1];
rz(-0.98140812) q[1];
sx q[1];
rz(2.9973105) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43898615) q[0];
sx q[0];
rz(-0.49285889) q[0];
sx q[0];
rz(2.4903557) q[0];
rz(-pi) q[1];
rz(2.2789196) q[2];
sx q[2];
rz(-1.7194334) q[2];
sx q[2];
rz(-1.3512063) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.49715091) q[1];
sx q[1];
rz(-1.9195917) q[1];
sx q[1];
rz(-1.3494874) q[1];
x q[2];
rz(-2.1142247) q[3];
sx q[3];
rz(-2.1680834) q[3];
sx q[3];
rz(-1.3940108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0353126) q[2];
sx q[2];
rz(-1.9929746) q[2];
sx q[2];
rz(-0.8379575) q[2];
rz(-2.641814) q[3];
sx q[3];
rz(-0.79334799) q[3];
sx q[3];
rz(0.70118538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(1.7734897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73117646) q[0];
sx q[0];
rz(-2.3623657) q[0];
sx q[0];
rz(0.91996361) q[0];
rz(-2.8255535) q[2];
sx q[2];
rz(-1.8506941) q[2];
sx q[2];
rz(-0.47376493) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4779799) q[1];
sx q[1];
rz(-0.77778947) q[1];
sx q[1];
rz(-1.6915093) q[1];
rz(2.1447324) q[3];
sx q[3];
rz(-1.878345) q[3];
sx q[3];
rz(2.779863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3086243) q[2];
sx q[2];
rz(-2.7699351) q[2];
sx q[2];
rz(-3.1404245) q[2];
rz(-0.33231169) q[3];
sx q[3];
rz(-1.6806335) q[3];
sx q[3];
rz(2.6619679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.371827) q[0];
sx q[0];
rz(-1.2170987) q[0];
sx q[0];
rz(2.0455072) q[0];
rz(-1.3990043) q[1];
sx q[1];
rz(-1.636248) q[1];
sx q[1];
rz(-0.91986626) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8815589) q[0];
sx q[0];
rz(-2.008201) q[0];
sx q[0];
rz(-3.0946671) q[0];
rz(-pi) q[1];
rz(-0.28025744) q[2];
sx q[2];
rz(-1.9748903) q[2];
sx q[2];
rz(1.0755444) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80650027) q[1];
sx q[1];
rz(-2.2119009) q[1];
sx q[1];
rz(-1.758746) q[1];
rz(-1.8125954) q[3];
sx q[3];
rz(-0.54587549) q[3];
sx q[3];
rz(2.161497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9225191) q[2];
sx q[2];
rz(-0.24524958) q[2];
sx q[2];
rz(1.5101439) q[2];
rz(-2.505693) q[3];
sx q[3];
rz(-0.70521516) q[3];
sx q[3];
rz(-2.6270134) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9566327) q[0];
sx q[0];
rz(-0.867046) q[0];
sx q[0];
rz(-1.4956723) q[0];
rz(3.0236859) q[1];
sx q[1];
rz(-1.3943744) q[1];
sx q[1];
rz(2.8848021) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.804931) q[0];
sx q[0];
rz(-1.1380702) q[0];
sx q[0];
rz(-3.0800729) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2725955) q[2];
sx q[2];
rz(-1.0062982) q[2];
sx q[2];
rz(0.80882458) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8578882) q[1];
sx q[1];
rz(-1.8562177) q[1];
sx q[1];
rz(-1.2976451) q[1];
x q[2];
rz(1.2076693) q[3];
sx q[3];
rz(-1.052891) q[3];
sx q[3];
rz(-2.7486441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7439277) q[2];
sx q[2];
rz(-2.3981018) q[2];
sx q[2];
rz(2.1916154) q[2];
rz(2.5395565) q[3];
sx q[3];
rz(-1.9046015) q[3];
sx q[3];
rz(-2.8034347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6947185) q[0];
sx q[0];
rz(-0.61114408) q[0];
sx q[0];
rz(-1.1993988) q[0];
rz(1.0097722) q[1];
sx q[1];
rz(-0.92242454) q[1];
sx q[1];
rz(-0.19933272) q[1];
rz(0.80712423) q[2];
sx q[2];
rz(-0.41584284) q[2];
sx q[2];
rz(-0.39932233) q[2];
rz(1.1789523) q[3];
sx q[3];
rz(-0.83189356) q[3];
sx q[3];
rz(-2.8883285) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
