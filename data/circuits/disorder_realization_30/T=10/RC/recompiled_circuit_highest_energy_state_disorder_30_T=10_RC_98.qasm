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
rz(0.93267814) q[0];
sx q[0];
rz(-2.7490766) q[0];
sx q[0];
rz(2.0445332) q[0];
rz(-1.8707844) q[1];
sx q[1];
rz(4.3391736) q[1];
sx q[1];
rz(14.270562) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2559833) q[0];
sx q[0];
rz(-0.75228158) q[0];
sx q[0];
rz(2.6153569) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8407048) q[2];
sx q[2];
rz(-0.2425293) q[2];
sx q[2];
rz(-2.1721942) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7157876) q[1];
sx q[1];
rz(-0.97570786) q[1];
sx q[1];
rz(-2.7462105) q[1];
rz(-pi) q[2];
rz(-0.013707073) q[3];
sx q[3];
rz(-0.64798149) q[3];
sx q[3];
rz(-1.5272087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2854332) q[2];
sx q[2];
rz(-2.9069052) q[2];
sx q[2];
rz(2.6402546) q[2];
rz(-1.8173789) q[3];
sx q[3];
rz(-0.54558498) q[3];
sx q[3];
rz(1.4668303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45601869) q[0];
sx q[0];
rz(-0.40559232) q[0];
sx q[0];
rz(1.4612041) q[0];
rz(-2.9623518) q[1];
sx q[1];
rz(-2.3724809) q[1];
sx q[1];
rz(-0.63757149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1792099) q[0];
sx q[0];
rz(-1.80733) q[0];
sx q[0];
rz(-1.3358227) q[0];
x q[1];
rz(-1.4316971) q[2];
sx q[2];
rz(-0.31965986) q[2];
sx q[2];
rz(-0.57397288) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5104622) q[1];
sx q[1];
rz(-2.3925214) q[1];
sx q[1];
rz(-2.1046348) q[1];
rz(-pi) q[2];
rz(-1.1944951) q[3];
sx q[3];
rz(-2.199389) q[3];
sx q[3];
rz(0.95910536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8403988) q[2];
sx q[2];
rz(-1.0777148) q[2];
sx q[2];
rz(3.0974498) q[2];
rz(-0.61331493) q[3];
sx q[3];
rz(-1.8200579) q[3];
sx q[3];
rz(0.80673748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.089207) q[0];
sx q[0];
rz(-3.0969924) q[0];
sx q[0];
rz(-1.3432304) q[0];
rz(1.7228458) q[1];
sx q[1];
rz(-2.2371465) q[1];
sx q[1];
rz(-0.80352965) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.270954) q[0];
sx q[0];
rz(-0.89193908) q[0];
sx q[0];
rz(-1.646433) q[0];
x q[1];
rz(-0.64068303) q[2];
sx q[2];
rz(-0.68924687) q[2];
sx q[2];
rz(-2.8169963) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.61345873) q[1];
sx q[1];
rz(-0.29859224) q[1];
sx q[1];
rz(-1.1996083) q[1];
rz(-pi) q[2];
rz(-0.66763287) q[3];
sx q[3];
rz(-2.5270695) q[3];
sx q[3];
rz(-0.4947084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.558305) q[2];
sx q[2];
rz(-2.1896157) q[2];
sx q[2];
rz(2.1236911) q[2];
rz(1.2244276) q[3];
sx q[3];
rz(-1.8095576) q[3];
sx q[3];
rz(-2.0126191) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1430436) q[0];
sx q[0];
rz(-2.3426549) q[0];
sx q[0];
rz(-2.1813188) q[0];
rz(-2.9472561) q[1];
sx q[1];
rz(-1.7002218) q[1];
sx q[1];
rz(0.0050553102) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95347217) q[0];
sx q[0];
rz(-1.0417099) q[0];
sx q[0];
rz(0.75624589) q[0];
rz(-pi) q[1];
rz(2.9580367) q[2];
sx q[2];
rz(-1.608299) q[2];
sx q[2];
rz(1.3958193) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9217164) q[1];
sx q[1];
rz(-1.7895234) q[1];
sx q[1];
rz(1.686386) q[1];
rz(-2.9937135) q[3];
sx q[3];
rz(-2.085146) q[3];
sx q[3];
rz(1.0837931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2527577) q[2];
sx q[2];
rz(-1.5450666) q[2];
sx q[2];
rz(-1.6952391) q[2];
rz(-1.8330005) q[3];
sx q[3];
rz(-1.7596217) q[3];
sx q[3];
rz(0.60477177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62110353) q[0];
sx q[0];
rz(-2.4172754) q[0];
sx q[0];
rz(-2.5526168) q[0];
rz(0.0079872459) q[1];
sx q[1];
rz(-2.0365448) q[1];
sx q[1];
rz(-2.0569107) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5863335) q[0];
sx q[0];
rz(-1.3166835) q[0];
sx q[0];
rz(1.3293367) q[0];
rz(-pi) q[1];
rz(0.27954526) q[2];
sx q[2];
rz(-1.5790006) q[2];
sx q[2];
rz(-1.4943124) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24655786) q[1];
sx q[1];
rz(-2.2854574) q[1];
sx q[1];
rz(-1.0393875) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6570377) q[3];
sx q[3];
rz(-0.6644333) q[3];
sx q[3];
rz(-2.2905615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7320554) q[2];
sx q[2];
rz(-1.0169225) q[2];
sx q[2];
rz(-1.3746877) q[2];
rz(0.095666766) q[3];
sx q[3];
rz(-1.5547662) q[3];
sx q[3];
rz(-2.5110551) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0365527) q[0];
sx q[0];
rz(-2.6281272) q[0];
sx q[0];
rz(1.8010944) q[0];
rz(2.4758677) q[1];
sx q[1];
rz(-1.3434854) q[1];
sx q[1];
rz(2.417477) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63987565) q[0];
sx q[0];
rz(-2.6611009) q[0];
sx q[0];
rz(2.2697422) q[0];
rz(-pi) q[1];
rz(-2.3409136) q[2];
sx q[2];
rz(-0.73298798) q[2];
sx q[2];
rz(-1.4633601) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.648702) q[1];
sx q[1];
rz(-0.48149432) q[1];
sx q[1];
rz(-2.8445811) q[1];
rz(3.136803) q[3];
sx q[3];
rz(-1.6671204) q[3];
sx q[3];
rz(-0.19118689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2035227) q[2];
sx q[2];
rz(-1.4504434) q[2];
sx q[2];
rz(-3.0435437) q[2];
rz(-0.34052643) q[3];
sx q[3];
rz(-2.2398658) q[3];
sx q[3];
rz(0.19671973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53511867) q[0];
sx q[0];
rz(-0.70347324) q[0];
sx q[0];
rz(3.086049) q[0];
rz(-0.37560383) q[1];
sx q[1];
rz(-0.77508488) q[1];
sx q[1];
rz(1.4449878) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7168769) q[0];
sx q[0];
rz(-2.3622241) q[0];
sx q[0];
rz(0.15044852) q[0];
rz(1.3169921) q[2];
sx q[2];
rz(-1.5198759) q[2];
sx q[2];
rz(1.6737991) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.58126175) q[1];
sx q[1];
rz(-0.37788299) q[1];
sx q[1];
rz(0.037200971) q[1];
x q[2];
rz(2.6747144) q[3];
sx q[3];
rz(-1.9345571) q[3];
sx q[3];
rz(-0.44565163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0982509) q[2];
sx q[2];
rz(-1.3037553) q[2];
sx q[2];
rz(-1.8275758) q[2];
rz(0.029953778) q[3];
sx q[3];
rz(-1.6311389) q[3];
sx q[3];
rz(2.8225074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7444262) q[0];
sx q[0];
rz(-0.63969669) q[0];
sx q[0];
rz(1.8901012) q[0];
rz(2.331612) q[1];
sx q[1];
rz(-0.73695838) q[1];
sx q[1];
rz(1.4553778) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43622323) q[0];
sx q[0];
rz(-1.3443265) q[0];
sx q[0];
rz(0.75948494) q[0];
rz(1.5157484) q[2];
sx q[2];
rz(-1.0491129) q[2];
sx q[2];
rz(-1.5417527) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4659459) q[1];
sx q[1];
rz(-2.1928804) q[1];
sx q[1];
rz(0.81333604) q[1];
rz(-pi) q[2];
rz(2.496893) q[3];
sx q[3];
rz(-0.94333269) q[3];
sx q[3];
rz(-2.6312466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.84888419) q[2];
sx q[2];
rz(-1.3564738) q[2];
sx q[2];
rz(1.9902309) q[2];
rz(0.014160841) q[3];
sx q[3];
rz(-1.3580946) q[3];
sx q[3];
rz(1.2667228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26361156) q[0];
sx q[0];
rz(-2.4284555) q[0];
sx q[0];
rz(2.1942595) q[0];
rz(0.5961279) q[1];
sx q[1];
rz(-1.7201741) q[1];
sx q[1];
rz(-1.5790342) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57347711) q[0];
sx q[0];
rz(-1.7422424) q[0];
sx q[0];
rz(-0.72400064) q[0];
rz(-pi) q[1];
x q[1];
rz(2.666108) q[2];
sx q[2];
rz(-1.4275996) q[2];
sx q[2];
rz(2.6729613) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4217675) q[1];
sx q[1];
rz(-0.9523069) q[1];
sx q[1];
rz(2.5773125) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72266717) q[3];
sx q[3];
rz(-2.8407466) q[3];
sx q[3];
rz(-2.2900555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8922213) q[2];
sx q[2];
rz(-0.5746848) q[2];
sx q[2];
rz(0.48386827) q[2];
rz(2.5441235) q[3];
sx q[3];
rz(-1.4474844) q[3];
sx q[3];
rz(-0.57815236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7396962) q[0];
sx q[0];
rz(-1.5902436) q[0];
sx q[0];
rz(0.36648146) q[0];
rz(-1.3260427) q[1];
sx q[1];
rz(-2.1791024) q[1];
sx q[1];
rz(2.0749626) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90049998) q[0];
sx q[0];
rz(-1.0783195) q[0];
sx q[0];
rz(-2.2482613) q[0];
rz(-pi) q[1];
rz(-2.5353955) q[2];
sx q[2];
rz(-1.6087039) q[2];
sx q[2];
rz(2.4647922) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6910256) q[1];
sx q[1];
rz(-0.98916173) q[1];
sx q[1];
rz(-2.8514991) q[1];
x q[2];
rz(-1.5068866) q[3];
sx q[3];
rz(-0.47807594) q[3];
sx q[3];
rz(1.687754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2610953) q[2];
sx q[2];
rz(-0.77793056) q[2];
sx q[2];
rz(2.0757389) q[2];
rz(-2.8737658) q[3];
sx q[3];
rz(-1.6949751) q[3];
sx q[3];
rz(2.6563787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5601226) q[0];
sx q[0];
rz(-2.1385834) q[0];
sx q[0];
rz(-2.9271097) q[0];
rz(-0.32196925) q[1];
sx q[1];
rz(-1.8245158) q[1];
sx q[1];
rz(1.239924) q[1];
rz(2.8602045) q[2];
sx q[2];
rz(-1.5783969) q[2];
sx q[2];
rz(1.6498298) q[2];
rz(-0.41245828) q[3];
sx q[3];
rz(-1.9720244) q[3];
sx q[3];
rz(0.26396863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
