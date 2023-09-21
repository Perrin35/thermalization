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
rz(-2.7789814) q[0];
rz(-2.2244722) q[1];
sx q[1];
rz(-2.6511104) q[1];
sx q[1];
rz(-2.7999556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65389079) q[0];
sx q[0];
rz(-0.39429769) q[0];
sx q[0];
rz(2.8173692) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50589675) q[2];
sx q[2];
rz(-0.6542754) q[2];
sx q[2];
rz(-1.8088532) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29400723) q[1];
sx q[1];
rz(-0.60677401) q[1];
sx q[1];
rz(1.0311544) q[1];
rz(-pi) q[2];
rz(0.2239286) q[3];
sx q[3];
rz(-1.180598) q[3];
sx q[3];
rz(-1.3003365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5477156) q[2];
sx q[2];
rz(-1.9888069) q[2];
sx q[2];
rz(0.95648742) q[2];
rz(-0.18125136) q[3];
sx q[3];
rz(-2.5528208) q[3];
sx q[3];
rz(1.6199002) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797121) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(-1.99615) q[0];
rz(-2.0727797) q[1];
sx q[1];
rz(-2.309598) q[1];
sx q[1];
rz(-0.72584814) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84859914) q[0];
sx q[0];
rz(-1.2147875) q[0];
sx q[0];
rz(-0.58649917) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5088596) q[2];
sx q[2];
rz(-0.69721141) q[2];
sx q[2];
rz(-2.7998507) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1314288) q[1];
sx q[1];
rz(-1.5457075) q[1];
sx q[1];
rz(1.8797092) q[1];
x q[2];
rz(-2.8511091) q[3];
sx q[3];
rz(-0.64290128) q[3];
sx q[3];
rz(0.089240616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6590865) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(2.144311) q[2];
rz(1.4128134) q[3];
sx q[3];
rz(-2.0480859) q[3];
sx q[3];
rz(2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41855758) q[0];
sx q[0];
rz(-0.96590531) q[0];
sx q[0];
rz(-2.0130656) q[0];
rz(1.8380802) q[1];
sx q[1];
rz(-2.4120657) q[1];
sx q[1];
rz(2.2361141) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826686) q[0];
sx q[0];
rz(-1.3805458) q[0];
sx q[0];
rz(0.24567901) q[0];
rz(-0.75240527) q[2];
sx q[2];
rz(-1.89626) q[2];
sx q[2];
rz(-2.4706555) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29618759) q[1];
sx q[1];
rz(-2.2803109) q[1];
sx q[1];
rz(-0.097464949) q[1];
x q[2];
rz(-0.043025322) q[3];
sx q[3];
rz(-0.85775162) q[3];
sx q[3];
rz(-0.42792861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7604312) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(-0.9872438) q[2];
rz(-3.0958214) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(-0.18019095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10953294) q[0];
sx q[0];
rz(-2.4592472) q[0];
sx q[0];
rz(-0.12606829) q[0];
rz(0.18445045) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(1.1674081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084598736) q[0];
sx q[0];
rz(-1.7718678) q[0];
sx q[0];
rz(-2.3720471) q[0];
x q[1];
rz(1.7502968) q[2];
sx q[2];
rz(-1.0738392) q[2];
sx q[2];
rz(-0.77979445) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3525225) q[1];
sx q[1];
rz(-1.8348502) q[1];
sx q[1];
rz(2.4697484) q[1];
x q[2];
rz(-2.3420742) q[3];
sx q[3];
rz(-0.82516731) q[3];
sx q[3];
rz(2.094401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2771153) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(2.5100822) q[2];
rz(-0.19550066) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(-0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9773848) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(0.50267977) q[0];
rz(1.1133105) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(2.6745093) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57468092) q[0];
sx q[0];
rz(-1.3710638) q[0];
sx q[0];
rz(1.1916222) q[0];
x q[1];
rz(-2.5555243) q[2];
sx q[2];
rz(-1.7826826) q[2];
sx q[2];
rz(0.52473247) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6131763) q[1];
sx q[1];
rz(-1.3532552) q[1];
sx q[1];
rz(-2.0030641) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2561982) q[3];
sx q[3];
rz(-1.6939031) q[3];
sx q[3];
rz(0.76243329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6749394) q[2];
sx q[2];
rz(-1.249524) q[2];
sx q[2];
rz(-0.49218991) q[2];
rz(-2.8594033) q[3];
sx q[3];
rz(-1.4331093) q[3];
sx q[3];
rz(-1.5658437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2714587) q[0];
sx q[0];
rz(-2.6014355) q[0];
sx q[0];
rz(0.085993275) q[0];
rz(1.4280691) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(1.496398) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2841227) q[0];
sx q[0];
rz(-2.5792482) q[0];
sx q[0];
rz(-1.1968632) q[0];
x q[1];
rz(-2.1093844) q[2];
sx q[2];
rz(-1.3628054) q[2];
sx q[2];
rz(-1.5185192) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8138258) q[1];
sx q[1];
rz(-1.5470978) q[1];
sx q[1];
rz(-1.9253255) q[1];
x q[2];
rz(0.71293998) q[3];
sx q[3];
rz(-1.561165) q[3];
sx q[3];
rz(2.9603017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66203403) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(0.39624828) q[2];
rz(0.59404343) q[3];
sx q[3];
rz(-1.0915353) q[3];
sx q[3];
rz(-1.5007639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8783766) q[0];
sx q[0];
rz(-2.5686869) q[0];
sx q[0];
rz(-2.0492045) q[0];
rz(1.3573525) q[1];
sx q[1];
rz(-1.7204294) q[1];
sx q[1];
rz(0.011627442) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8788293) q[0];
sx q[0];
rz(-1.3063523) q[0];
sx q[0];
rz(2.7468365) q[0];
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
x q[0];
rz(-1.1639451) q[1];
sx q[1];
rz(-1.5944591) q[1];
sx q[1];
rz(2.1169099) q[1];
x q[2];
rz(-0.39832468) q[3];
sx q[3];
rz(-0.71361226) q[3];
sx q[3];
rz(-0.45809612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.45550436) q[2];
sx q[2];
rz(-2.4839165) q[2];
sx q[2];
rz(0.28798506) q[2];
rz(-1.1503495) q[3];
sx q[3];
rz(-1.8201613) q[3];
sx q[3];
rz(-2.5434125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31036723) q[0];
sx q[0];
rz(-0.53124017) q[0];
sx q[0];
rz(-1.9294552) q[0];
rz(-0.37777004) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(-1.4935965) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15570116) q[0];
sx q[0];
rz(-0.37030664) q[0];
sx q[0];
rz(2.5168688) q[0];
x q[1];
rz(3.0817501) q[2];
sx q[2];
rz(-1.4388196) q[2];
sx q[2];
rz(3.0795385) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89477506) q[1];
sx q[1];
rz(-0.90118876) q[1];
sx q[1];
rz(2.0991904) q[1];
rz(3.0231608) q[3];
sx q[3];
rz(-1.0716972) q[3];
sx q[3];
rz(-2.0407608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7887855) q[2];
sx q[2];
rz(-1.1151423) q[2];
sx q[2];
rz(1.8898233) q[2];
rz(-2.2552323) q[3];
sx q[3];
rz(-0.76615196) q[3];
sx q[3];
rz(3.0322976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733474) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(2.9633203) q[0];
rz(3.0687304) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(-1.1791621) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68547738) q[0];
sx q[0];
rz(-0.59887409) q[0];
sx q[0];
rz(-2.6045538) q[0];
rz(-pi) q[1];
rz(-0.91782848) q[2];
sx q[2];
rz(-0.73790109) q[2];
sx q[2];
rz(-0.81673056) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33680962) q[1];
sx q[1];
rz(-2.0840624) q[1];
sx q[1];
rz(1.8670765) q[1];
rz(-0.9411373) q[3];
sx q[3];
rz(-2.5459873) q[3];
sx q[3];
rz(-0.60048088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4420085) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(-1.0317831) q[2];
rz(1.7992841) q[3];
sx q[3];
rz(-1.7248025) q[3];
sx q[3];
rz(0.22527307) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28829065) q[0];
sx q[0];
rz(-2.0993711) q[0];
sx q[0];
rz(-0.061696079) q[0];
rz(-1.306698) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(1.0888938) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93931224) q[0];
sx q[0];
rz(-1.6334403) q[0];
sx q[0];
rz(1.5059727) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91572362) q[2];
sx q[2];
rz(-1.7172126) q[2];
sx q[2];
rz(2.5953948) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2606925) q[1];
sx q[1];
rz(-2.0339194) q[1];
sx q[1];
rz(-2.5093964) q[1];
x q[2];
rz(-1.0763361) q[3];
sx q[3];
rz(-2.7571602) q[3];
sx q[3];
rz(-1.9264551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5237727) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(-0.069996746) q[2];
rz(-2.2648515) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4927647) q[0];
sx q[0];
rz(-1.6236826) q[0];
sx q[0];
rz(0.59326011) q[0];
rz(1.6246673) q[1];
sx q[1];
rz(-1.0472052) q[1];
sx q[1];
rz(-0.44890961) q[1];
rz(-0.91296997) q[2];
sx q[2];
rz(-1.2428478) q[2];
sx q[2];
rz(-1.7075677) q[2];
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