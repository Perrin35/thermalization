OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1113623) q[0];
sx q[0];
rz(-2.4863939) q[0];
sx q[0];
rz(-0.97638786) q[0];
rz(-0.916565) q[1];
sx q[1];
rz(2.38148) q[1];
sx q[1];
rz(10.217477) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9206031) q[0];
sx q[0];
rz(-1.5467637) q[0];
sx q[0];
rz(-0.8997957) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7723299) q[2];
sx q[2];
rz(-2.5662176) q[2];
sx q[2];
rz(-1.0550261) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77289171) q[1];
sx q[1];
rz(-1.4979616) q[1];
sx q[1];
rz(-1.1198977) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3000852) q[3];
sx q[3];
rz(-2.2268647) q[3];
sx q[3];
rz(2.5049202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4684326) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(2.6426278) q[2];
rz(0.4237825) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(-3.1288778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38133165) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(0.99200845) q[0];
rz(-0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-2.7761249) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1680981) q[0];
sx q[0];
rz(-1.1143648) q[0];
sx q[0];
rz(-1.2903851) q[0];
rz(-pi) q[1];
rz(-2.2720049) q[2];
sx q[2];
rz(-0.80539942) q[2];
sx q[2];
rz(2.373901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9531627) q[1];
sx q[1];
rz(-2.1663482) q[1];
sx q[1];
rz(-1.1577391) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5378546) q[3];
sx q[3];
rz(-2.9085277) q[3];
sx q[3];
rz(0.92249289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3084597) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(-1.0773405) q[2];
rz(-2.9789553) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(-1.2619654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4804374) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(-2.5307632) q[0];
rz(1.707533) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(-2.6909713) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47942802) q[0];
sx q[0];
rz(-3.056884) q[0];
sx q[0];
rz(-1.9730294) q[0];
x q[1];
rz(0.30655105) q[2];
sx q[2];
rz(-1.9286641) q[2];
sx q[2];
rz(-1.4150261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0011117) q[1];
sx q[1];
rz(-1.2621659) q[1];
sx q[1];
rz(0.66092296) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87617989) q[3];
sx q[3];
rz(-1.3324705) q[3];
sx q[3];
rz(-0.5660457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11015636) q[2];
sx q[2];
rz(-1.0855731) q[2];
sx q[2];
rz(-1.9499367) q[2];
rz(0.84447652) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0464756) q[0];
sx q[0];
rz(-1.3866383) q[0];
sx q[0];
rz(-2.0003831) q[0];
rz(-1.2966688) q[1];
sx q[1];
rz(-2.7099687) q[1];
sx q[1];
rz(-0.30532125) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6396219) q[0];
sx q[0];
rz(-2.1272215) q[0];
sx q[0];
rz(1.200838) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4112169) q[2];
sx q[2];
rz(-2.5502898) q[2];
sx q[2];
rz(-3.0808466) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0381283) q[1];
sx q[1];
rz(-2.3410019) q[1];
sx q[1];
rz(-0.14726463) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3577945) q[3];
sx q[3];
rz(-2.4664306) q[3];
sx q[3];
rz(-2.6019707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9332463) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(-2.7944881) q[2];
rz(2.7582205) q[3];
sx q[3];
rz(-0.88630668) q[3];
sx q[3];
rz(-2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7326996) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(0.53891671) q[0];
rz(-1.4096889) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(0.40571037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.273542) q[0];
sx q[0];
rz(-1.8450071) q[0];
sx q[0];
rz(-2.4540841) q[0];
rz(1.1666127) q[2];
sx q[2];
rz(-2.5378072) q[2];
sx q[2];
rz(0.052978901) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87863641) q[1];
sx q[1];
rz(-2.7107781) q[1];
sx q[1];
rz(-2.2104435) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3375072) q[3];
sx q[3];
rz(-2.3217215) q[3];
sx q[3];
rz(-1.6398721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6695909) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(1.6873138) q[2];
rz(2.3029095) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(1.7680761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8307761) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(0.40711656) q[0];
rz(2.282417) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(-1.8242594) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2989527) q[0];
sx q[0];
rz(-1.5484122) q[0];
sx q[0];
rz(1.0277261) q[0];
x q[1];
rz(0.69465722) q[2];
sx q[2];
rz(-2.3431871) q[2];
sx q[2];
rz(0.77381221) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6588905) q[1];
sx q[1];
rz(-1.5735978) q[1];
sx q[1];
rz(0.83562327) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7139211) q[3];
sx q[3];
rz(-1.0276405) q[3];
sx q[3];
rz(0.021051858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.09981) q[2];
sx q[2];
rz(-2.2618497) q[2];
sx q[2];
rz(0.89522925) q[2];
rz(1.9334531) q[3];
sx q[3];
rz(-0.67966214) q[3];
sx q[3];
rz(0.7081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5773425) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(2.0136925) q[0];
rz(0.212542) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(-1.2316661) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82176498) q[0];
sx q[0];
rz(-1.2152831) q[0];
sx q[0];
rz(-0.54152821) q[0];
rz(-pi) q[1];
rz(-0.77881323) q[2];
sx q[2];
rz(-2.3328569) q[2];
sx q[2];
rz(0.29701172) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1612138) q[1];
sx q[1];
rz(-1.8002138) q[1];
sx q[1];
rz(1.6769665) q[1];
x q[2];
rz(-1.6024186) q[3];
sx q[3];
rz(-0.60472371) q[3];
sx q[3];
rz(-1.5550169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3125375) q[2];
sx q[2];
rz(-0.91784358) q[2];
sx q[2];
rz(-2.5404239) q[2];
rz(-0.51856315) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(1.7140088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.9940051) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(-3.0684379) q[0];
rz(-1.0271429) q[1];
sx q[1];
rz(-1.606753) q[1];
sx q[1];
rz(0.59246078) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1661467) q[0];
sx q[0];
rz(-1.6172292) q[0];
sx q[0];
rz(-2.5509044) q[0];
rz(-1.0560889) q[2];
sx q[2];
rz(-1.0564431) q[2];
sx q[2];
rz(-2.8714542) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29443024) q[1];
sx q[1];
rz(-1.8995598) q[1];
sx q[1];
rz(-0.86818236) q[1];
rz(-pi) q[2];
rz(-0.27505625) q[3];
sx q[3];
rz(-1.961879) q[3];
sx q[3];
rz(-1.0712136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.56117326) q[2];
sx q[2];
rz(-1.4658804) q[2];
sx q[2];
rz(0.44211659) q[2];
rz(2.2413065) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(0.013307868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(3.1135547) q[0];
sx q[0];
rz(-0.91911626) q[0];
sx q[0];
rz(-1.8369209) q[0];
rz(1.6944983) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(2.5794199) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32414831) q[0];
sx q[0];
rz(-1.2153017) q[0];
sx q[0];
rz(2.4118511) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.047692) q[2];
sx q[2];
rz(-2.0996089) q[2];
sx q[2];
rz(1.7319861) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8852946) q[1];
sx q[1];
rz(-2.2513299) q[1];
sx q[1];
rz(2.2018432) q[1];
rz(-0.77769827) q[3];
sx q[3];
rz(-2.5986528) q[3];
sx q[3];
rz(2.267595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8582981) q[2];
sx q[2];
rz(-1.7389745) q[2];
sx q[2];
rz(0.98158681) q[2];
rz(-3.1372519) q[3];
sx q[3];
rz(-1.0467632) q[3];
sx q[3];
rz(0.59345746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8269862) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(0.04709588) q[0];
rz(1.3866407) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(3.093739) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78794107) q[0];
sx q[0];
rz(-1.5997868) q[0];
sx q[0];
rz(-0.91584648) q[0];
rz(-pi) q[1];
rz(-1.4638897) q[2];
sx q[2];
rz(-0.44471834) q[2];
sx q[2];
rz(0.61483785) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.88356) q[1];
sx q[1];
rz(-2.0000534) q[1];
sx q[1];
rz(-2.2116823) q[1];
x q[2];
rz(-1.6258532) q[3];
sx q[3];
rz(-1.1686472) q[3];
sx q[3];
rz(-1.0226585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64426595) q[2];
sx q[2];
rz(-1.1844119) q[2];
sx q[2];
rz(-2.2052374) q[2];
rz(-0.80091536) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(2.463533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8235648) q[0];
sx q[0];
rz(-2.2021273) q[0];
sx q[0];
rz(0.20903023) q[0];
rz(-0.50280747) q[1];
sx q[1];
rz(-1.3192758) q[1];
sx q[1];
rz(-2.1154138) q[1];
rz(-1.5488831) q[2];
sx q[2];
rz(-1.0196601) q[2];
sx q[2];
rz(-1.7658726) q[2];
rz(-2.7754178) q[3];
sx q[3];
rz(-1.1911285) q[3];
sx q[3];
rz(1.7366684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
