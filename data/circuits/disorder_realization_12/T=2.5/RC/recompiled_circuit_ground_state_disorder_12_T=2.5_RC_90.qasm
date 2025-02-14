OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1801017) q[0];
sx q[0];
rz(-0.56668007) q[0];
sx q[0];
rz(-0.79071796) q[0];
rz(-0.92791954) q[1];
sx q[1];
rz(-0.035049573) q[1];
sx q[1];
rz(-2.8383281) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7105616) q[0];
sx q[0];
rz(-0.91791422) q[0];
sx q[0];
rz(1.2827355) q[0];
rz(-pi) q[1];
rz(0.82178564) q[2];
sx q[2];
rz(-0.50210929) q[2];
sx q[2];
rz(-2.6977476) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3825079) q[1];
sx q[1];
rz(-2.0469991) q[1];
sx q[1];
rz(0.97530968) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6536775) q[3];
sx q[3];
rz(-2.304416) q[3];
sx q[3];
rz(0.26700156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7683893) q[2];
sx q[2];
rz(-2.781227) q[2];
sx q[2];
rz(-2.4599794) q[2];
rz(-2.5810589) q[3];
sx q[3];
rz(-2.3571641) q[3];
sx q[3];
rz(0.64422977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19175567) q[0];
sx q[0];
rz(-2.2332709) q[0];
sx q[0];
rz(0.35582304) q[0];
rz(-0.25335723) q[1];
sx q[1];
rz(-2.2919877) q[1];
sx q[1];
rz(-1.8960948) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6424774) q[0];
sx q[0];
rz(-1.0986975) q[0];
sx q[0];
rz(1.8224707) q[0];
rz(-pi) q[1];
rz(-0.63388692) q[2];
sx q[2];
rz(-1.1340464) q[2];
sx q[2];
rz(2.2417807) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51735866) q[1];
sx q[1];
rz(-1.9962203) q[1];
sx q[1];
rz(0.58601554) q[1];
rz(-pi) q[2];
rz(0.50470501) q[3];
sx q[3];
rz(-1.0796121) q[3];
sx q[3];
rz(-2.1757656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7443098) q[2];
sx q[2];
rz(-0.77703589) q[2];
sx q[2];
rz(3.0139319) q[2];
rz(-0.67576659) q[3];
sx q[3];
rz(-0.45857576) q[3];
sx q[3];
rz(2.7202386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0449988) q[0];
sx q[0];
rz(-0.24462044) q[0];
sx q[0];
rz(2.007572) q[0];
rz(0.89316142) q[1];
sx q[1];
rz(-0.90407073) q[1];
sx q[1];
rz(1.8697416) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24090996) q[0];
sx q[0];
rz(-2.5782452) q[0];
sx q[0];
rz(-1.6202089) q[0];
rz(-pi) q[1];
x q[1];
rz(1.261146) q[2];
sx q[2];
rz(-1.4952588) q[2];
sx q[2];
rz(0.47020082) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.73619702) q[1];
sx q[1];
rz(-0.73976019) q[1];
sx q[1];
rz(-1.2642047) q[1];
x q[2];
rz(-0.16678129) q[3];
sx q[3];
rz(-1.3540917) q[3];
sx q[3];
rz(-3.0971017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21358718) q[2];
sx q[2];
rz(-2.2760133) q[2];
sx q[2];
rz(1.340284) q[2];
rz(-1.853893) q[3];
sx q[3];
rz(-1.4475334) q[3];
sx q[3];
rz(0.48937669) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61557788) q[0];
sx q[0];
rz(-2.5647793) q[0];
sx q[0];
rz(2.4082129) q[0];
rz(-2.3461657) q[1];
sx q[1];
rz(-0.19618244) q[1];
sx q[1];
rz(-1.9733852) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0709658) q[0];
sx q[0];
rz(-2.7160108) q[0];
sx q[0];
rz(1.312567) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57046095) q[2];
sx q[2];
rz(-2.1242363) q[2];
sx q[2];
rz(-1.2177856) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.296098) q[1];
sx q[1];
rz(-1.5278399) q[1];
sx q[1];
rz(1.6011493) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9549035) q[3];
sx q[3];
rz(-2.1400937) q[3];
sx q[3];
rz(-0.74654859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9613793) q[2];
sx q[2];
rz(-2.3258379) q[2];
sx q[2];
rz(-0.95480314) q[2];
rz(2.09156) q[3];
sx q[3];
rz(-2.3539216) q[3];
sx q[3];
rz(2.5900456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10057218) q[0];
sx q[0];
rz(-2.6028778) q[0];
sx q[0];
rz(2.4464497) q[0];
rz(1.9255385) q[1];
sx q[1];
rz(-1.0168409) q[1];
sx q[1];
rz(0.020811828) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3580768) q[0];
sx q[0];
rz(-2.4672271) q[0];
sx q[0];
rz(1.3410617) q[0];
rz(-pi) q[1];
rz(-0.32692744) q[2];
sx q[2];
rz(-2.280664) q[2];
sx q[2];
rz(2.8500789) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4592308) q[1];
sx q[1];
rz(-1.8098847) q[1];
sx q[1];
rz(-3.1296631) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0486879) q[3];
sx q[3];
rz(-1.6321401) q[3];
sx q[3];
rz(-2.7829704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9713328) q[2];
sx q[2];
rz(-2.8016165) q[2];
sx q[2];
rz(-0.60410947) q[2];
rz(-0.64770925) q[3];
sx q[3];
rz(-2.6205687) q[3];
sx q[3];
rz(2.6807396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0172417) q[0];
sx q[0];
rz(-1.0931953) q[0];
sx q[0];
rz(-2.7830615) q[0];
rz(2.5784946) q[1];
sx q[1];
rz(-0.9022572) q[1];
sx q[1];
rz(-2.8360352) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53154778) q[0];
sx q[0];
rz(-2.1912337) q[0];
sx q[0];
rz(-2.0187747) q[0];
x q[1];
rz(1.3995869) q[2];
sx q[2];
rz(-1.7504901) q[2];
sx q[2];
rz(3.0269738) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7737732) q[1];
sx q[1];
rz(-1.0350643) q[1];
sx q[1];
rz(2.1526835) q[1];
x q[2];
rz(-0.077353954) q[3];
sx q[3];
rz(-1.6408444) q[3];
sx q[3];
rz(-1.9590247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2629338) q[2];
sx q[2];
rz(-0.96076751) q[2];
sx q[2];
rz(0.052138694) q[2];
rz(0.42929286) q[3];
sx q[3];
rz(-1.7310127) q[3];
sx q[3];
rz(-0.2521635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4572064) q[0];
sx q[0];
rz(-2.6173213) q[0];
sx q[0];
rz(-2.9065409) q[0];
rz(-2.5382407) q[1];
sx q[1];
rz(-2.3124606) q[1];
sx q[1];
rz(-2.5650909) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2430778) q[0];
sx q[0];
rz(-0.91204005) q[0];
sx q[0];
rz(-2.2284449) q[0];
rz(-pi) q[1];
rz(-2.7965661) q[2];
sx q[2];
rz(-0.80586159) q[2];
sx q[2];
rz(-2.8553183) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53945069) q[1];
sx q[1];
rz(-1.2509402) q[1];
sx q[1];
rz(1.1509239) q[1];
rz(-1.0888388) q[3];
sx q[3];
rz(-1.4432943) q[3];
sx q[3];
rz(-1.3393928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0779886) q[2];
sx q[2];
rz(-0.84090191) q[2];
sx q[2];
rz(-1.3464751) q[2];
rz(0.50955647) q[3];
sx q[3];
rz(-1.2600803) q[3];
sx q[3];
rz(-0.30663651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52617306) q[0];
sx q[0];
rz(-2.058448) q[0];
sx q[0];
rz(0.35588595) q[0];
rz(0.7016167) q[1];
sx q[1];
rz(-2.9569148) q[1];
sx q[1];
rz(0.96431771) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0052781) q[0];
sx q[0];
rz(-1.1169254) q[0];
sx q[0];
rz(-0.14880609) q[0];
x q[1];
rz(-0.020267296) q[2];
sx q[2];
rz(-1.8635443) q[2];
sx q[2];
rz(-2.0890638) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0605392) q[1];
sx q[1];
rz(-2.6463228) q[1];
sx q[1];
rz(3.1002863) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4889433) q[3];
sx q[3];
rz(-1.8813895) q[3];
sx q[3];
rz(1.5378584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.648916) q[2];
sx q[2];
rz(-0.70768386) q[2];
sx q[2];
rz(-0.81563449) q[2];
rz(-2.2260769) q[3];
sx q[3];
rz(-0.86796498) q[3];
sx q[3];
rz(-0.22667949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8409519) q[0];
sx q[0];
rz(-2.9673567) q[0];
sx q[0];
rz(1.2409644) q[0];
rz(0.10494431) q[1];
sx q[1];
rz(-0.34067708) q[1];
sx q[1];
rz(-2.6822283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37840415) q[0];
sx q[0];
rz(-1.3306466) q[0];
sx q[0];
rz(1.2021007) q[0];
x q[1];
rz(2.8723889) q[2];
sx q[2];
rz(-1.891231) q[2];
sx q[2];
rz(1.073369) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8018559) q[1];
sx q[1];
rz(-1.1366399) q[1];
sx q[1];
rz(-1.3375907) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9271803) q[3];
sx q[3];
rz(-2.5567856) q[3];
sx q[3];
rz(2.6406276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0910054) q[2];
sx q[2];
rz(-2.4531187) q[2];
sx q[2];
rz(-3.0509994) q[2];
rz(-0.76720864) q[3];
sx q[3];
rz(-1.6588666) q[3];
sx q[3];
rz(1.4513133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7200658) q[0];
sx q[0];
rz(-0.92362112) q[0];
sx q[0];
rz(-1.5371171) q[0];
rz(-2.1859258) q[1];
sx q[1];
rz(-1.4612528) q[1];
sx q[1];
rz(-0.54735267) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5794649) q[0];
sx q[0];
rz(-1.5648995) q[0];
sx q[0];
rz(1.6505169) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59835206) q[2];
sx q[2];
rz(-0.6171591) q[2];
sx q[2];
rz(-1.1197108) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0136614) q[1];
sx q[1];
rz(-1.8056925) q[1];
sx q[1];
rz(-2.3460745) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0738952) q[3];
sx q[3];
rz(-0.92313952) q[3];
sx q[3];
rz(-1.5941509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1110111) q[2];
sx q[2];
rz(-0.17685282) q[2];
sx q[2];
rz(2.492823) q[2];
rz(-2.9099921) q[3];
sx q[3];
rz(-0.88445556) q[3];
sx q[3];
rz(3.054936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(-1.9301013) q[0];
sx q[0];
rz(1.7420266) q[0];
rz(0.65921417) q[1];
sx q[1];
rz(-1.2792239) q[1];
sx q[1];
rz(-1.4469133) q[1];
rz(1.7991015) q[2];
sx q[2];
rz(-1.6189499) q[2];
sx q[2];
rz(1.5789643) q[2];
rz(0.054604758) q[3];
sx q[3];
rz(-1.7163897) q[3];
sx q[3];
rz(-0.39672273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
