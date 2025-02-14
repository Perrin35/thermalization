OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73622048) q[0];
sx q[0];
rz(3.1625746) q[0];
sx q[0];
rz(7.4641135) q[0];
rz(2.6612072) q[1];
sx q[1];
rz(-1.1552224) q[1];
sx q[1];
rz(0.46673271) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8202756) q[0];
sx q[0];
rz(-1.2877687) q[0];
sx q[0];
rz(-0.77518605) q[0];
rz(-pi) q[1];
rz(2.4619815) q[2];
sx q[2];
rz(-2.4386897) q[2];
sx q[2];
rz(1.694569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8267104) q[1];
sx q[1];
rz(-2.9032384) q[1];
sx q[1];
rz(1.5895859) q[1];
x q[2];
rz(-1.0084413) q[3];
sx q[3];
rz(-2.1690024) q[3];
sx q[3];
rz(1.4149208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0537609) q[2];
sx q[2];
rz(-1.6873282) q[2];
sx q[2];
rz(-1.9161179) q[2];
rz(-0.37442225) q[3];
sx q[3];
rz(-1.1332847) q[3];
sx q[3];
rz(-2.5882914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18749172) q[0];
sx q[0];
rz(-0.88045374) q[0];
sx q[0];
rz(-0.53400293) q[0];
rz(-0.39132896) q[1];
sx q[1];
rz(-0.44328872) q[1];
sx q[1];
rz(2.2543529) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92838642) q[0];
sx q[0];
rz(-1.22661) q[0];
sx q[0];
rz(1.2437938) q[0];
x q[1];
rz(0.75669719) q[2];
sx q[2];
rz(-1.9644418) q[2];
sx q[2];
rz(1.4243038) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5767867) q[1];
sx q[1];
rz(-2.8999355) q[1];
sx q[1];
rz(2.5741626) q[1];
rz(2.7999185) q[3];
sx q[3];
rz(-1.0782584) q[3];
sx q[3];
rz(1.1600174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5947764) q[2];
sx q[2];
rz(-1.3204601) q[2];
sx q[2];
rz(-2.7776264) q[2];
rz(-1.4173896) q[3];
sx q[3];
rz(-0.28984362) q[3];
sx q[3];
rz(-2.3072306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2367547) q[0];
sx q[0];
rz(-0.098243864) q[0];
sx q[0];
rz(-2.8114317) q[0];
rz(-0.58369058) q[1];
sx q[1];
rz(-2.3287562) q[1];
sx q[1];
rz(2.6469753) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9512105) q[0];
sx q[0];
rz(-1.3000164) q[0];
sx q[0];
rz(-0.27373224) q[0];
rz(-pi) q[1];
rz(1.0801267) q[2];
sx q[2];
rz(-1.3441836) q[2];
sx q[2];
rz(-1.8267711) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3266625) q[1];
sx q[1];
rz(-2.0771793) q[1];
sx q[1];
rz(-0.5270919) q[1];
x q[2];
rz(0.98683896) q[3];
sx q[3];
rz(-1.47746) q[3];
sx q[3];
rz(-0.10658857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3459449) q[2];
sx q[2];
rz(-2.4730885) q[2];
sx q[2];
rz(-0.12269679) q[2];
rz(-3.0986541) q[3];
sx q[3];
rz(-2.0624845) q[3];
sx q[3];
rz(0.19449657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4075459) q[0];
sx q[0];
rz(-2.6341944) q[0];
sx q[0];
rz(-0.88974446) q[0];
rz(-0.45097688) q[1];
sx q[1];
rz(-0.86800066) q[1];
sx q[1];
rz(2.6873592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4047336) q[0];
sx q[0];
rz(-1.6772207) q[0];
sx q[0];
rz(2.7989796) q[0];
rz(2.7142482) q[2];
sx q[2];
rz(-1.5691535) q[2];
sx q[2];
rz(-2.9438007) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52773038) q[1];
sx q[1];
rz(-2.3017028) q[1];
sx q[1];
rz(1.1682603) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20619456) q[3];
sx q[3];
rz(-0.36199328) q[3];
sx q[3];
rz(0.21807204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.98895994) q[2];
sx q[2];
rz(-0.6344513) q[2];
sx q[2];
rz(-2.3839942) q[2];
rz(0.18975137) q[3];
sx q[3];
rz(-1.8228143) q[3];
sx q[3];
rz(-0.54722133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2106638) q[0];
sx q[0];
rz(-2.3543816) q[0];
sx q[0];
rz(1.8481365) q[0];
rz(-2.5594607) q[1];
sx q[1];
rz(-1.4385834) q[1];
sx q[1];
rz(2.9621946) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2820601) q[0];
sx q[0];
rz(-0.27320394) q[0];
sx q[0];
rz(0.84635107) q[0];
rz(-pi) q[1];
rz(1.8151692) q[2];
sx q[2];
rz(-2.2269911) q[2];
sx q[2];
rz(0.49733053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.071956858) q[1];
sx q[1];
rz(-1.7577802) q[1];
sx q[1];
rz(2.2025639) q[1];
x q[2];
rz(0.50197451) q[3];
sx q[3];
rz(-0.93794146) q[3];
sx q[3];
rz(-2.8090734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77357117) q[2];
sx q[2];
rz(-2.7715235) q[2];
sx q[2];
rz(0.5698815) q[2];
rz(-2.701345) q[3];
sx q[3];
rz(-1.8972242) q[3];
sx q[3];
rz(-2.7888035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2920947) q[0];
sx q[0];
rz(-2.8773913) q[0];
sx q[0];
rz(2.0000892) q[0];
rz(-1.3453206) q[1];
sx q[1];
rz(-0.99015403) q[1];
sx q[1];
rz(0.17123953) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2889338) q[0];
sx q[0];
rz(-1.6717413) q[0];
sx q[0];
rz(-2.6600718) q[0];
rz(0.87821853) q[2];
sx q[2];
rz(-0.64413749) q[2];
sx q[2];
rz(-0.53976166) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0192193) q[1];
sx q[1];
rz(-1.7369441) q[1];
sx q[1];
rz(-2.5801587) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1275131) q[3];
sx q[3];
rz(-1.5601592) q[3];
sx q[3];
rz(-2.6680846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3693927) q[2];
sx q[2];
rz(-0.93595305) q[2];
sx q[2];
rz(-3.0041223) q[2];
rz(1.2482268) q[3];
sx q[3];
rz(-0.56709254) q[3];
sx q[3];
rz(-1.9198157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(0.084994706) q[0];
sx q[0];
rz(-2.5039112) q[0];
sx q[0];
rz(1.6531264) q[0];
rz(0.78530637) q[1];
sx q[1];
rz(-1.7410802) q[1];
sx q[1];
rz(-1.8280425) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5564329) q[0];
sx q[0];
rz(-1.4476416) q[0];
sx q[0];
rz(0.086104546) q[0];
x q[1];
rz(1.4573757) q[2];
sx q[2];
rz(-1.7608425) q[2];
sx q[2];
rz(-2.6148877) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0848799) q[1];
sx q[1];
rz(-0.093357714) q[1];
sx q[1];
rz(2.8356524) q[1];
rz(-pi) q[2];
rz(1.63229) q[3];
sx q[3];
rz(-2.4447421) q[3];
sx q[3];
rz(1.9887672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5588348) q[2];
sx q[2];
rz(-1.02966) q[2];
sx q[2];
rz(2.6666857) q[2];
rz(2.9389935) q[3];
sx q[3];
rz(-2.30862) q[3];
sx q[3];
rz(-1.1420265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.1982034) q[0];
sx q[0];
rz(-3.0084963) q[0];
sx q[0];
rz(-0.30580795) q[0];
rz(0.43931475) q[1];
sx q[1];
rz(-2.3190934) q[1];
sx q[1];
rz(1.3154715) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5852642) q[0];
sx q[0];
rz(-1.5595115) q[0];
sx q[0];
rz(0.68345921) q[0];
rz(-0.28355174) q[2];
sx q[2];
rz(-2.1075948) q[2];
sx q[2];
rz(-2.5187522) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.55637344) q[1];
sx q[1];
rz(-1.6917768) q[1];
sx q[1];
rz(0.78775265) q[1];
rz(-2.1602116) q[3];
sx q[3];
rz(-1.9187315) q[3];
sx q[3];
rz(-0.98435005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7266015) q[2];
sx q[2];
rz(-1.5982268) q[2];
sx q[2];
rz(0.42465633) q[2];
rz(-3.1096544) q[3];
sx q[3];
rz(-1.5141124) q[3];
sx q[3];
rz(-0.64938515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53751078) q[0];
sx q[0];
rz(-0.6518971) q[0];
sx q[0];
rz(-0.55484581) q[0];
rz(2.8765053) q[1];
sx q[1];
rz(-1.4684497) q[1];
sx q[1];
rz(1.2010942) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5145668) q[0];
sx q[0];
rz(-1.1794588) q[0];
sx q[0];
rz(-2.5364801) q[0];
rz(-0.55055203) q[2];
sx q[2];
rz(-0.32062396) q[2];
sx q[2];
rz(-0.52832802) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9399413) q[1];
sx q[1];
rz(-1.9186416) q[1];
sx q[1];
rz(-1.6099206) q[1];
rz(-pi) q[2];
rz(0.74984925) q[3];
sx q[3];
rz(-1.5886652) q[3];
sx q[3];
rz(-0.32040473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.25192866) q[2];
sx q[2];
rz(-0.98968518) q[2];
sx q[2];
rz(0.84452334) q[2];
rz(1.1318413) q[3];
sx q[3];
rz(-0.90513888) q[3];
sx q[3];
rz(-1.3593146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.1691386) q[0];
sx q[0];
rz(-0.3903946) q[0];
sx q[0];
rz(-1.9688695) q[0];
rz(-2.4523465) q[1];
sx q[1];
rz(-2.0447562) q[1];
sx q[1];
rz(-2.8776339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64484064) q[0];
sx q[0];
rz(-2.1979316) q[0];
sx q[0];
rz(-1.0154614) q[0];
rz(-0.9969292) q[2];
sx q[2];
rz(-2.3341114) q[2];
sx q[2];
rz(1.6107744) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.86706272) q[1];
sx q[1];
rz(-1.2559051) q[1];
sx q[1];
rz(-2.9462183) q[1];
rz(0.12216025) q[3];
sx q[3];
rz(-0.93046832) q[3];
sx q[3];
rz(-1.2608918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3802203) q[2];
sx q[2];
rz(-1.994588) q[2];
sx q[2];
rz(2.0210361) q[2];
rz(-1.606696) q[3];
sx q[3];
rz(-2.1963019) q[3];
sx q[3];
rz(-1.5103316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071038889) q[0];
sx q[0];
rz(-2.4865535) q[0];
sx q[0];
rz(3.0336663) q[0];
rz(0.72036605) q[1];
sx q[1];
rz(-1.9279059) q[1];
sx q[1];
rz(2.7005213) q[1];
rz(0.82577827) q[2];
sx q[2];
rz(-2.0407531) q[2];
sx q[2];
rz(-2.4498418) q[2];
rz(0.10590774) q[3];
sx q[3];
rz(-0.61951588) q[3];
sx q[3];
rz(-2.9168152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
