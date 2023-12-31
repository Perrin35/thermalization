OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5126098) q[0];
sx q[0];
rz(-1.0273758) q[0];
sx q[0];
rz(0.36261121) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(2.6511104) q[1];
sx q[1];
rz(9.766415) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1385961) q[0];
sx q[0];
rz(-1.9435104) q[0];
sx q[0];
rz(1.4390104) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59092893) q[2];
sx q[2];
rz(-1.2714296) q[2];
sx q[2];
rz(0.17609827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8475854) q[1];
sx q[1];
rz(-2.5348186) q[1];
sx q[1];
rz(-1.0311544) q[1];
x q[2];
rz(-1.0757252) q[3];
sx q[3];
rz(-2.6945811) q[3];
sx q[3];
rz(-1.8398374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5477156) q[2];
sx q[2];
rz(-1.9888069) q[2];
sx q[2];
rz(2.1851052) q[2];
rz(-0.18125136) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(1.5216924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797121) q[0];
sx q[0];
rz(-3.0759838) q[0];
sx q[0];
rz(1.99615) q[0];
rz(1.068813) q[1];
sx q[1];
rz(-2.309598) q[1];
sx q[1];
rz(-0.72584814) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84859914) q[0];
sx q[0];
rz(-1.9268052) q[0];
sx q[0];
rz(2.5550935) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63273301) q[2];
sx q[2];
rz(-0.69721141) q[2];
sx q[2];
rz(2.7998507) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.010163807) q[1];
sx q[1];
rz(-1.5457075) q[1];
sx q[1];
rz(-1.8797092) q[1];
x q[2];
rz(-2.5190982) q[3];
sx q[3];
rz(-1.3982292) q[3];
sx q[3];
rz(-1.7164001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6590865) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(-2.144311) q[2];
rz(-1.7287792) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.7230351) q[0];
sx q[0];
rz(-2.1756873) q[0];
sx q[0];
rz(-1.1285271) q[0];
rz(1.8380802) q[1];
sx q[1];
rz(-2.4120657) q[1];
sx q[1];
rz(-0.9054786) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3826686) q[0];
sx q[0];
rz(-1.3805458) q[0];
sx q[0];
rz(2.8959136) q[0];
x q[1];
rz(2.682914) q[2];
sx q[2];
rz(-2.3346666) q[2];
sx q[2];
rz(-0.57070953) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.14721522) q[1];
sx q[1];
rz(-2.4265687) q[1];
sx q[1];
rz(-1.6836402) q[1];
rz(0.85729349) q[3];
sx q[3];
rz(-1.5382574) q[3];
sx q[3];
rz(-1.1710222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3811615) q[2];
sx q[2];
rz(-2.9734549) q[2];
sx q[2];
rz(0.9872438) q[2];
rz(-0.045771249) q[3];
sx q[3];
rz(-1.8493303) q[3];
sx q[3];
rz(-0.18019095) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0320597) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(-0.12606829) q[0];
rz(-2.9571422) q[1];
sx q[1];
rz(-1.0686864) q[1];
sx q[1];
rz(1.9741845) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0569939) q[0];
sx q[0];
rz(-1.3697249) q[0];
sx q[0];
rz(-2.3720471) q[0];
rz(-pi) q[1];
rz(1.3912958) q[2];
sx q[2];
rz(-1.0738392) q[2];
sx q[2];
rz(-2.3617982) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.78907012) q[1];
sx q[1];
rz(-1.3067424) q[1];
sx q[1];
rz(-0.67184429) q[1];
x q[2];
rz(0.91058369) q[3];
sx q[3];
rz(-1.0331717) q[3];
sx q[3];
rz(1.1066574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8644774) q[2];
sx q[2];
rz(-0.41431674) q[2];
sx q[2];
rz(-0.63151044) q[2];
rz(-0.19550066) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(-0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
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
rz(-1.9773848) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(-0.50267977) q[0];
rz(1.1133105) q[1];
sx q[1];
rz(-1.0031507) q[1];
sx q[1];
rz(-2.6745093) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0665849) q[0];
sx q[0];
rz(-1.1995312) q[0];
sx q[0];
rz(-2.9270372) q[0];
rz(-2.7706625) q[2];
sx q[2];
rz(-2.5226454) q[2];
sx q[2];
rz(-1.3528454) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.621532) q[1];
sx q[1];
rz(-2.660775) q[1];
sx q[1];
rz(2.0562999) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3777556) q[3];
sx q[3];
rz(-2.4469864) q[3];
sx q[3];
rz(0.65929268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.46665329) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(2.6494027) q[2];
rz(-0.28218937) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(-1.5658437) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2714587) q[0];
sx q[0];
rz(-0.54015714) q[0];
sx q[0];
rz(3.0555994) q[0];
rz(1.7135235) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(-1.496398) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7487885) q[0];
sx q[0];
rz(-1.7668056) q[0];
sx q[0];
rz(-1.0402354) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9005269) q[2];
sx q[2];
rz(-1.0450372) q[2];
sx q[2];
rz(0.17503967) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3277668) q[1];
sx q[1];
rz(-1.5944949) q[1];
sx q[1];
rz(-1.9253255) q[1];
rz(3.1268678) q[3];
sx q[3];
rz(-2.428599) q[3];
sx q[3];
rz(-1.7409489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66203403) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(2.7453444) q[2];
rz(0.59404343) q[3];
sx q[3];
rz(-1.0915353) q[3];
sx q[3];
rz(1.6408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2632161) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(1.0923882) q[0];
rz(1.3573525) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(3.1299652) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8893338) q[0];
sx q[0];
rz(-0.47124915) q[0];
sx q[0];
rz(0.61347368) q[0];
x q[1];
rz(-1.5989223) q[2];
sx q[2];
rz(-1.5453494) q[2];
sx q[2];
rz(0.59288073) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7736518) q[1];
sx q[1];
rz(-2.5950187) q[1];
sx q[1];
rz(-1.6163338) q[1];
x q[2];
rz(1.8947951) q[3];
sx q[3];
rz(-2.218459) q[3];
sx q[3];
rz(-0.049829359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.45550436) q[2];
sx q[2];
rz(-2.4839165) q[2];
sx q[2];
rz(2.8536076) q[2];
rz(1.1503495) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(-2.5434125) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31036723) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(1.9294552) q[0];
rz(-2.7638226) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(-1.6479962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.134682) q[0];
sx q[0];
rz(-1.784076) q[0];
sx q[0];
rz(-2.8365305) q[0];
rz(-pi) q[1];
rz(1.703007) q[2];
sx q[2];
rz(-1.5114748) q[2];
sx q[2];
rz(1.5166264) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8131866) q[1];
sx q[1];
rz(-1.9771736) q[1];
sx q[1];
rz(-2.3996668) q[1];
x q[2];
rz(-1.0687374) q[3];
sx q[3];
rz(-1.6747253) q[3];
sx q[3];
rz(-2.7285189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7887855) q[2];
sx q[2];
rz(-1.1151423) q[2];
sx q[2];
rz(1.8898233) q[2];
rz(-0.88636032) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(3.0322976) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66824526) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(-2.9633203) q[0];
rz(-0.072862236) q[1];
sx q[1];
rz(-0.59385121) q[1];
sx q[1];
rz(1.1791621) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68547738) q[0];
sx q[0];
rz(-2.5427186) q[0];
sx q[0];
rz(2.6045538) q[0];
x q[1];
rz(0.5046919) q[2];
sx q[2];
rz(-2.1345277) q[2];
sx q[2];
rz(-1.6187514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33680962) q[1];
sx q[1];
rz(-2.0840624) q[1];
sx q[1];
rz(-1.2745162) q[1];
rz(-pi) q[2];
rz(2.0719028) q[3];
sx q[3];
rz(-1.9074829) q[3];
sx q[3];
rz(1.5130373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4420085) q[2];
sx q[2];
rz(-2.1506491) q[2];
sx q[2];
rz(1.0317831) q[2];
rz(-1.3423086) q[3];
sx q[3];
rz(-1.7248025) q[3];
sx q[3];
rz(-2.9163196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28829065) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(-3.0798966) q[0];
rz(1.8348947) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(-2.0526989) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2022804) q[0];
sx q[0];
rz(-1.5081524) q[0];
sx q[0];
rz(-1.63562) q[0];
rz(-pi) q[1];
rz(-0.18386545) q[2];
sx q[2];
rz(-2.2176748) q[2];
sx q[2];
rz(-2.0053787) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1354462) q[1];
sx q[1];
rz(-1.0137614) q[1];
sx q[1];
rz(-2.1250493) q[1];
rz(2.0652566) q[3];
sx q[3];
rz(-0.38443243) q[3];
sx q[3];
rz(1.9264551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6178199) q[2];
sx q[2];
rz(-0.68796316) q[2];
sx q[2];
rz(-3.0715959) q[2];
rz(-0.87674117) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(-2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6488279) q[0];
sx q[0];
rz(-1.5179101) q[0];
sx q[0];
rz(-2.5483325) q[0];
rz(-1.5169253) q[1];
sx q[1];
rz(-1.0472052) q[1];
sx q[1];
rz(-0.44890961) q[1];
rz(2.7355315) q[2];
sx q[2];
rz(-0.95352298) q[2];
sx q[2];
rz(2.7609115) q[2];
rz(-1.2354479) q[3];
sx q[3];
rz(-1.609758) q[3];
sx q[3];
rz(2.5555425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
