OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6289829) q[0];
sx q[0];
rz(-2.1142168) q[0];
sx q[0];
rz(2.7789814) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(-0.4904823) q[1];
sx q[1];
rz(-0.34163707) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65389079) q[0];
sx q[0];
rz(-0.39429769) q[0];
sx q[0];
rz(-0.32422347) q[0];
rz(-pi) q[1];
rz(-0.59092893) q[2];
sx q[2];
rz(-1.2714296) q[2];
sx q[2];
rz(-2.9654944) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.33597782) q[1];
sx q[1];
rz(-1.059638) q[1];
sx q[1];
rz(-2.7989945) q[1];
rz(0.2239286) q[3];
sx q[3];
rz(-1.180598) q[3];
sx q[3];
rz(-1.3003365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5477156) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(-0.95648742) q[2];
rz(2.9603413) q[3];
sx q[3];
rz(-2.5528208) q[3];
sx q[3];
rz(-1.5216924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0618806) q[0];
sx q[0];
rz(-3.0759838) q[0];
sx q[0];
rz(-1.99615) q[0];
rz(-2.0727797) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(-2.4157445) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2929935) q[0];
sx q[0];
rz(-1.9268052) q[0];
sx q[0];
rz(2.5550935) q[0];
rz(-pi) q[1];
rz(2.5088596) q[2];
sx q[2];
rz(-2.4443812) q[2];
sx q[2];
rz(2.7998507) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4821856) q[1];
sx q[1];
rz(-2.8316951) q[1];
sx q[1];
rz(1.6531497) q[1];
x q[2];
rz(-2.8511091) q[3];
sx q[3];
rz(-2.4986914) q[3];
sx q[3];
rz(3.052352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6590865) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(2.144311) q[2];
rz(1.7287792) q[3];
sx q[3];
rz(-1.0935067) q[3];
sx q[3];
rz(2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(-1.8380802) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(2.2361141) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826686) q[0];
sx q[0];
rz(-1.3805458) q[0];
sx q[0];
rz(-2.8959136) q[0];
rz(-2.682914) q[2];
sx q[2];
rz(-2.3346666) q[2];
sx q[2];
rz(-2.5708831) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.29618759) q[1];
sx q[1];
rz(-2.2803109) q[1];
sx q[1];
rz(0.097464949) q[1];
x q[2];
rz(-2.2842992) q[3];
sx q[3];
rz(-1.6033353) q[3];
sx q[3];
rz(-1.9705704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.7604312) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(0.9872438) q[2];
rz(0.045771249) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(2.9614017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10953294) q[0];
sx q[0];
rz(-2.4592472) q[0];
sx q[0];
rz(3.0155244) q[0];
rz(-2.9571422) q[1];
sx q[1];
rz(-1.0686864) q[1];
sx q[1];
rz(1.9741845) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0569939) q[0];
sx q[0];
rz(-1.3697249) q[0];
sx q[0];
rz(0.76954557) q[0];
x q[1];
rz(-2.6378176) q[2];
sx q[2];
rz(-1.7283895) q[2];
sx q[2];
rz(0.7047082) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5645204) q[1];
sx q[1];
rz(-2.2153691) q[1];
sx q[1];
rz(1.2381899) q[1];
rz(-pi) q[2];
rz(2.3420742) q[3];
sx q[3];
rz(-0.82516731) q[3];
sx q[3];
rz(1.0471917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2771153) q[2];
sx q[2];
rz(-0.41431674) q[2];
sx q[2];
rz(-2.5100822) q[2];
rz(0.19550066) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9773848) q[0];
sx q[0];
rz(-0.13810869) q[0];
sx q[0];
rz(-0.50267977) q[0];
rz(-2.0282822) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(2.6745093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5669117) q[0];
sx q[0];
rz(-1.7705288) q[0];
sx q[0];
rz(1.1916222) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5555243) q[2];
sx q[2];
rz(-1.7826826) q[2];
sx q[2];
rz(-2.6168602) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6131763) q[1];
sx q[1];
rz(-1.3532552) q[1];
sx q[1];
rz(1.1385285) q[1];
rz(-pi) q[2];
rz(-1.3777556) q[3];
sx q[3];
rz(-0.69460624) q[3];
sx q[3];
rz(2.4823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46665329) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(2.6494027) q[2];
rz(0.28218937) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(-1.5757489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2714587) q[0];
sx q[0];
rz(-2.6014355) q[0];
sx q[0];
rz(0.085993275) q[0];
rz(1.4280691) q[1];
sx q[1];
rz(-2.1336887) q[1];
sx q[1];
rz(1.6451947) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85747) q[0];
sx q[0];
rz(-2.5792482) q[0];
sx q[0];
rz(1.1968632) q[0];
rz(-2.1093844) q[2];
sx q[2];
rz(-1.7787873) q[2];
sx q[2];
rz(-1.6230735) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8138258) q[1];
sx q[1];
rz(-1.5470978) q[1];
sx q[1];
rz(1.9253255) q[1];
rz(-2.4286527) q[3];
sx q[3];
rz(-1.561165) q[3];
sx q[3];
rz(2.9603017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.66203403) q[2];
sx q[2];
rz(-2.5223314) q[2];
sx q[2];
rz(0.39624828) q[2];
rz(-2.5475492) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(-1.6408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783766) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(-2.0492045) q[0];
rz(-1.7842402) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(3.1299652) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7250925) q[0];
sx q[0];
rz(-1.19048) q[0];
sx q[0];
rz(-1.8561383) q[0];
rz(-pi) q[1];
rz(-0.83519148) q[2];
sx q[2];
rz(-3.1036658) q[2];
sx q[2];
rz(1.4284301) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7736518) q[1];
sx q[1];
rz(-0.546574) q[1];
sx q[1];
rz(1.5252588) q[1];
rz(-pi) q[2];
rz(0.67354789) q[3];
sx q[3];
rz(-1.3140956) q[3];
sx q[3];
rz(1.7208769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6860883) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(2.8536076) q[2];
rz(1.1503495) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(0.59818017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31036723) q[0];
sx q[0];
rz(-0.53124017) q[0];
sx q[0];
rz(1.2121375) q[0];
rz(0.37777004) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(1.4935965) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6389248) q[0];
sx q[0];
rz(-1.272861) q[0];
sx q[0];
rz(-1.7940679) q[0];
x q[1];
rz(1.703007) q[2];
sx q[2];
rz(-1.5114748) q[2];
sx q[2];
rz(1.5166264) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.89477506) q[1];
sx q[1];
rz(-2.2404039) q[1];
sx q[1];
rz(-1.0424022) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7842403) q[3];
sx q[3];
rz(-0.5118013) q[3];
sx q[3];
rz(-1.7970999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7887855) q[2];
sx q[2];
rz(-2.0264503) q[2];
sx q[2];
rz(1.8898233) q[2];
rz(-2.2552323) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(-3.0322976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66824526) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(0.1782724) q[0];
rz(-0.072862236) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(-1.1791621) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0807063) q[0];
sx q[0];
rz(-2.076425) q[0];
sx q[0];
rz(-1.9067184) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2237642) q[2];
sx q[2];
rz(-0.73790109) q[2];
sx q[2];
rz(0.81673056) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.804783) q[1];
sx q[1];
rz(-1.0575302) q[1];
sx q[1];
rz(-1.8670765) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7618802) q[3];
sx q[3];
rz(-2.0413997) q[3];
sx q[3];
rz(3.0203366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4420085) q[2];
sx q[2];
rz(-2.1506491) q[2];
sx q[2];
rz(2.1098095) q[2];
rz(1.3423086) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(-2.9163196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.853302) q[0];
sx q[0];
rz(-2.0993711) q[0];
sx q[0];
rz(3.0798966) q[0];
rz(-1.306698) q[1];
sx q[1];
rz(-1.4790165) q[1];
sx q[1];
rz(2.0526989) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93931224) q[0];
sx q[0];
rz(-1.6334403) q[0];
sx q[0];
rz(-1.5059727) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3332974) q[2];
sx q[2];
rz(-2.4727159) q[2];
sx q[2];
rz(0.83692688) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.006146487) q[1];
sx q[1];
rz(-2.1278312) q[1];
sx q[1];
rz(2.1250493) q[1];
x q[2];
rz(1.2286931) q[3];
sx q[3];
rz(-1.3918687) q[3];
sx q[3];
rz(0.10781328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6178199) q[2];
sx q[2];
rz(-0.68796316) q[2];
sx q[2];
rz(-3.0715959) q[2];
rz(-2.2648515) q[3];
sx q[3];
rz(-1.3249967) q[3];
sx q[3];
rz(0.35818067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.6488279) q[0];
sx q[0];
rz(-1.5179101) q[0];
sx q[0];
rz(-2.5483325) q[0];
rz(1.5169253) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(1.0629874) q[2];
sx q[2];
rz(-2.4175736) q[2];
sx q[2];
rz(-2.8833817) q[2];
rz(0.041257507) q[3];
sx q[3];
rz(-1.2357124) q[3];
sx q[3];
rz(0.99832051) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];