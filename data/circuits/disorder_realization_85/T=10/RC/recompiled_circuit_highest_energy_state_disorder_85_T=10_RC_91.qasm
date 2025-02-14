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
rz(0.99211168) q[0];
sx q[0];
rz(3.8881128) q[0];
sx q[0];
rz(9.9915656) q[0];
rz(1.2504638) q[1];
sx q[1];
rz(-1.992978) q[1];
sx q[1];
rz(2.221938) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5772668) q[0];
sx q[0];
rz(-0.79123215) q[0];
sx q[0];
rz(1.7971695) q[0];
x q[1];
rz(-2.4162499) q[2];
sx q[2];
rz(-1.1769466) q[2];
sx q[2];
rz(2.3462636) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2055605) q[1];
sx q[1];
rz(-1.524393) q[1];
sx q[1];
rz(-0.64179365) q[1];
rz(-pi) q[2];
rz(-1.7149431) q[3];
sx q[3];
rz(-2.1021772) q[3];
sx q[3];
rz(-1.2856158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.93996843) q[2];
sx q[2];
rz(-2.2087966) q[2];
sx q[2];
rz(-2.0882108) q[2];
rz(-2.4273704) q[3];
sx q[3];
rz(-2.768399) q[3];
sx q[3];
rz(1.8478954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5193704) q[0];
sx q[0];
rz(-2.433233) q[0];
sx q[0];
rz(-2.569662) q[0];
rz(-1.8427303) q[1];
sx q[1];
rz(-0.79890257) q[1];
sx q[1];
rz(-2.3518708) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2970726) q[0];
sx q[0];
rz(-1.4747341) q[0];
sx q[0];
rz(1.8734832) q[0];
x q[1];
rz(2.4873494) q[2];
sx q[2];
rz(-1.6983319) q[2];
sx q[2];
rz(2.5728284) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4724628) q[1];
sx q[1];
rz(-1.8505368) q[1];
sx q[1];
rz(0.6196687) q[1];
rz(-0.67482194) q[3];
sx q[3];
rz(-1.1742419) q[3];
sx q[3];
rz(0.048585437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.85084891) q[2];
sx q[2];
rz(-0.76883832) q[2];
sx q[2];
rz(1.7791465) q[2];
rz(2.8507774) q[3];
sx q[3];
rz(-1.7751866) q[3];
sx q[3];
rz(-0.10663685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0027851) q[0];
sx q[0];
rz(-1.0951575) q[0];
sx q[0];
rz(2.441067) q[0];
rz(-2.6558212) q[1];
sx q[1];
rz(-2.3120717) q[1];
sx q[1];
rz(0.22451678) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6922263) q[0];
sx q[0];
rz(-2.7590519) q[0];
sx q[0];
rz(-0.40479779) q[0];
x q[1];
rz(-0.098578886) q[2];
sx q[2];
rz(-0.99462992) q[2];
sx q[2];
rz(2.2599932) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3923334) q[1];
sx q[1];
rz(-0.9888923) q[1];
sx q[1];
rz(0.090163377) q[1];
x q[2];
rz(2.2994479) q[3];
sx q[3];
rz(-1.4858467) q[3];
sx q[3];
rz(2.2865975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38480467) q[2];
sx q[2];
rz(-0.89185682) q[2];
sx q[2];
rz(-3.0925114) q[2];
rz(-0.067042025) q[3];
sx q[3];
rz(-1.1578355) q[3];
sx q[3];
rz(-2.4750347) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9662358) q[0];
sx q[0];
rz(-0.55713621) q[0];
sx q[0];
rz(0.019158451) q[0];
rz(-1.7823904) q[1];
sx q[1];
rz(-1.4298341) q[1];
sx q[1];
rz(-0.0044936831) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86165262) q[0];
sx q[0];
rz(-0.77095448) q[0];
sx q[0];
rz(2.3907135) q[0];
x q[1];
rz(-2.4158813) q[2];
sx q[2];
rz(-2.5134676) q[2];
sx q[2];
rz(-0.84838644) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7178065) q[1];
sx q[1];
rz(-0.39925925) q[1];
sx q[1];
rz(-2.0171793) q[1];
rz(-pi) q[2];
rz(1.6636623) q[3];
sx q[3];
rz(-1.7735274) q[3];
sx q[3];
rz(-2.7331238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4517453) q[2];
sx q[2];
rz(-1.0901901) q[2];
sx q[2];
rz(-1.6881662) q[2];
rz(1.2545741) q[3];
sx q[3];
rz(-1.6202319) q[3];
sx q[3];
rz(-0.5622676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8101863) q[0];
sx q[0];
rz(-1.9047381) q[0];
sx q[0];
rz(-1.9796665) q[0];
rz(-2.3792073) q[1];
sx q[1];
rz(-0.98680174) q[1];
sx q[1];
rz(-0.42957482) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9305161) q[0];
sx q[0];
rz(-1.1628224) q[0];
sx q[0];
rz(0.039647722) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75865443) q[2];
sx q[2];
rz(-1.7650371) q[2];
sx q[2];
rz(1.5264896) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.030368806) q[1];
sx q[1];
rz(-2.3850523) q[1];
sx q[1];
rz(-1.0683505) q[1];
x q[2];
rz(0.91752538) q[3];
sx q[3];
rz(-1.8442196) q[3];
sx q[3];
rz(1.7961674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21542159) q[2];
sx q[2];
rz(-1.2913707) q[2];
sx q[2];
rz(1.8298836) q[2];
rz(2.6808776) q[3];
sx q[3];
rz(-1.0034794) q[3];
sx q[3];
rz(-3.0084012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8120414) q[0];
sx q[0];
rz(-2.9819745) q[0];
sx q[0];
rz(2.0352236) q[0];
rz(0.1144935) q[1];
sx q[1];
rz(-2.0888927) q[1];
sx q[1];
rz(0.053650275) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44941266) q[0];
sx q[0];
rz(-2.8859038) q[0];
sx q[0];
rz(-1.5935672) q[0];
rz(-pi) q[1];
rz(-1.6752376) q[2];
sx q[2];
rz(-0.76600961) q[2];
sx q[2];
rz(-3.0867689) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1785894) q[1];
sx q[1];
rz(-0.51437639) q[1];
sx q[1];
rz(-1.6806755) q[1];
rz(-pi) q[2];
rz(2.1838004) q[3];
sx q[3];
rz(-1.0292064) q[3];
sx q[3];
rz(-1.047618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8657118) q[2];
sx q[2];
rz(-1.2578473) q[2];
sx q[2];
rz(2.6344521) q[2];
rz(-1.7670613) q[3];
sx q[3];
rz(-1.5406698) q[3];
sx q[3];
rz(1.1486294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62992612) q[0];
sx q[0];
rz(-2.0138795) q[0];
sx q[0];
rz(-0.041286904) q[0];
rz(-2.4033026) q[1];
sx q[1];
rz(-1.9169044) q[1];
sx q[1];
rz(2.1521177) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7558003) q[0];
sx q[0];
rz(-0.72052252) q[0];
sx q[0];
rz(2.3260443) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7858753) q[2];
sx q[2];
rz(-1.748744) q[2];
sx q[2];
rz(-1.2426113) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7478024) q[1];
sx q[1];
rz(-0.89266864) q[1];
sx q[1];
rz(-2.3403779) q[1];
rz(-pi) q[2];
rz(2.2253898) q[3];
sx q[3];
rz(-1.8713017) q[3];
sx q[3];
rz(2.253807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5993293) q[2];
sx q[2];
rz(-2.0928536) q[2];
sx q[2];
rz(-0.85912022) q[2];
rz(3.1298992) q[3];
sx q[3];
rz(-1.499736) q[3];
sx q[3];
rz(-0.11277994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3149253) q[0];
sx q[0];
rz(-2.5496917) q[0];
sx q[0];
rz(-2.9097606) q[0];
rz(0.58206093) q[1];
sx q[1];
rz(-1.8128017) q[1];
sx q[1];
rz(-1.9225165) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7340379) q[0];
sx q[0];
rz(-2.1269607) q[0];
sx q[0];
rz(1.4329877) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6232331) q[2];
sx q[2];
rz(-1.2182256) q[2];
sx q[2];
rz(0.8964311) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.2208913) q[1];
sx q[1];
rz(-0.24020837) q[1];
sx q[1];
rz(0.88418737) q[1];
rz(-pi) q[2];
rz(2.3858504) q[3];
sx q[3];
rz(-1.983194) q[3];
sx q[3];
rz(2.2568767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9570738) q[2];
sx q[2];
rz(-1.3513214) q[2];
sx q[2];
rz(-0.78682023) q[2];
rz(1.6400853) q[3];
sx q[3];
rz(-2.7101176) q[3];
sx q[3];
rz(-1.4260346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.96104923) q[0];
sx q[0];
rz(-1.3732055) q[0];
sx q[0];
rz(-0.94026646) q[0];
rz(2.6967948) q[1];
sx q[1];
rz(-0.85591379) q[1];
sx q[1];
rz(-2.99517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0068374182) q[0];
sx q[0];
rz(-1.869258) q[0];
sx q[0];
rz(-1.7025856) q[0];
x q[1];
rz(1.8249885) q[2];
sx q[2];
rz(-1.8349525) q[2];
sx q[2];
rz(-0.025042695) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.9616016) q[1];
sx q[1];
rz(-2.4117984) q[1];
sx q[1];
rz(-0.37288937) q[1];
rz(-1.4692471) q[3];
sx q[3];
rz(-2.3595105) q[3];
sx q[3];
rz(2.8080432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.08708295) q[2];
sx q[2];
rz(-1.6528218) q[2];
sx q[2];
rz(-0.83756891) q[2];
rz(2.2086823) q[3];
sx q[3];
rz(-0.98740238) q[3];
sx q[3];
rz(-2.3605409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26302108) q[0];
sx q[0];
rz(-1.2282547) q[0];
sx q[0];
rz(1.3680869) q[0];
rz(3.0655762) q[1];
sx q[1];
rz(-2.5612505) q[1];
sx q[1];
rz(-2.0215633) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90864414) q[0];
sx q[0];
rz(-0.41260037) q[0];
sx q[0];
rz(2.2792675) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70856382) q[2];
sx q[2];
rz(-1.3363133) q[2];
sx q[2];
rz(-1.3402779) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8260086) q[1];
sx q[1];
rz(-0.30696973) q[1];
sx q[1];
rz(1.2493285) q[1];
x q[2];
rz(0.10346966) q[3];
sx q[3];
rz(-1.7333687) q[3];
sx q[3];
rz(-1.181447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7398305) q[2];
sx q[2];
rz(-1.6900475) q[2];
sx q[2];
rz(2.8945727) q[2];
rz(-0.57010993) q[3];
sx q[3];
rz(-2.6542122) q[3];
sx q[3];
rz(1.1183687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0982672) q[0];
sx q[0];
rz(-1.7541616) q[0];
sx q[0];
rz(2.7761205) q[0];
rz(-0.098516057) q[1];
sx q[1];
rz(-2.5010074) q[1];
sx q[1];
rz(0.88190257) q[1];
rz(-1.1578887) q[2];
sx q[2];
rz(-1.5831686) q[2];
sx q[2];
rz(-1.8148212) q[2];
rz(-2.7612272) q[3];
sx q[3];
rz(-1.813253) q[3];
sx q[3];
rz(-1.062275) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
