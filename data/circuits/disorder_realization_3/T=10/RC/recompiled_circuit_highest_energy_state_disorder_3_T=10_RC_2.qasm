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
rz(-2.6920707) q[0];
sx q[0];
rz(-1.7188526) q[0];
sx q[0];
rz(-1.0499522) q[0];
rz(2.6226251) q[1];
sx q[1];
rz(-1.608404) q[1];
sx q[1];
rz(-0.050921507) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22171709) q[0];
sx q[0];
rz(-1.1394986) q[0];
sx q[0];
rz(-1.9274516) q[0];
rz(-pi) q[1];
rz(1.2363288) q[2];
sx q[2];
rz(-1.7750472) q[2];
sx q[2];
rz(-1.1676015) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9999034) q[1];
sx q[1];
rz(-1.2517271) q[1];
sx q[1];
rz(1.2525108) q[1];
rz(-1.8597569) q[3];
sx q[3];
rz(-0.70631344) q[3];
sx q[3];
rz(3.0970124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.27429399) q[2];
sx q[2];
rz(-2.2653502) q[2];
sx q[2];
rz(1.743861) q[2];
rz(2.6869669) q[3];
sx q[3];
rz(-2.2635098) q[3];
sx q[3];
rz(-1.4286058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6468663) q[0];
sx q[0];
rz(-0.88748256) q[0];
sx q[0];
rz(-0.60321641) q[0];
rz(1.8869205) q[1];
sx q[1];
rz(-2.5277977) q[1];
sx q[1];
rz(-0.57993728) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3361084) q[0];
sx q[0];
rz(-1.0382129) q[0];
sx q[0];
rz(1.7642154) q[0];
rz(-pi) q[1];
rz(-2.0236778) q[2];
sx q[2];
rz(-1.0476026) q[2];
sx q[2];
rz(2.4459237) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3160657) q[1];
sx q[1];
rz(-2.5361119) q[1];
sx q[1];
rz(1.5010631) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5744561) q[3];
sx q[3];
rz(-1.8923762) q[3];
sx q[3];
rz(2.0547607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4846399) q[2];
sx q[2];
rz(-2.3886949) q[2];
sx q[2];
rz(0.37772712) q[2];
rz(-1.7083302) q[3];
sx q[3];
rz(-1.5092311) q[3];
sx q[3];
rz(1.9506955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3483873) q[0];
sx q[0];
rz(-3.0866525) q[0];
sx q[0];
rz(-1.6435664) q[0];
rz(-1.9801697) q[1];
sx q[1];
rz(-1.252251) q[1];
sx q[1];
rz(2.2983671) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4068309) q[0];
sx q[0];
rz(-2.1929682) q[0];
sx q[0];
rz(0.0733331) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21768985) q[2];
sx q[2];
rz(-1.161631) q[2];
sx q[2];
rz(0.37693757) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23318849) q[1];
sx q[1];
rz(-0.68906765) q[1];
sx q[1];
rz(-1.3702964) q[1];
rz(-pi) q[2];
rz(2.4405582) q[3];
sx q[3];
rz(-1.8314519) q[3];
sx q[3];
rz(3.1150027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2735542) q[2];
sx q[2];
rz(-0.66323438) q[2];
sx q[2];
rz(0.794945) q[2];
rz(0.76977175) q[3];
sx q[3];
rz(-2.218518) q[3];
sx q[3];
rz(-2.9845089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.58761) q[0];
sx q[0];
rz(-1.3251323) q[0];
sx q[0];
rz(-0.28833589) q[0];
rz(-0.68710697) q[1];
sx q[1];
rz(-1.4930875) q[1];
sx q[1];
rz(-1.4289325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9014578) q[0];
sx q[0];
rz(-0.013738886) q[0];
sx q[0];
rz(0.49679784) q[0];
rz(0.32055118) q[2];
sx q[2];
rz(-0.42141576) q[2];
sx q[2];
rz(1.6089862) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.3521139) q[1];
sx q[1];
rz(-1.5055471) q[1];
sx q[1];
rz(-2.8852374) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5023584) q[3];
sx q[3];
rz(-2.1393288) q[3];
sx q[3];
rz(-0.45445874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.571542) q[2];
sx q[2];
rz(-0.96633458) q[2];
sx q[2];
rz(0.16560444) q[2];
rz(-0.064519493) q[3];
sx q[3];
rz(-0.12972984) q[3];
sx q[3];
rz(1.8701514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9128543) q[0];
sx q[0];
rz(-1.9485291) q[0];
sx q[0];
rz(0.91113973) q[0];
rz(0.97081026) q[1];
sx q[1];
rz(-1.8152922) q[1];
sx q[1];
rz(-0.86311805) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2902888) q[0];
sx q[0];
rz(-0.57803854) q[0];
sx q[0];
rz(2.1208288) q[0];
rz(-pi) q[1];
rz(-1.6550894) q[2];
sx q[2];
rz(-0.26940036) q[2];
sx q[2];
rz(1.0625372) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.94416241) q[1];
sx q[1];
rz(-1.3384377) q[1];
sx q[1];
rz(0.54411035) q[1];
x q[2];
rz(-3.0961995) q[3];
sx q[3];
rz(-1.8828585) q[3];
sx q[3];
rz(1.6873716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0866278) q[2];
sx q[2];
rz(-1.8885771) q[2];
sx q[2];
rz(-0.48913726) q[2];
rz(2.1431811) q[3];
sx q[3];
rz(-0.17397927) q[3];
sx q[3];
rz(-0.099099549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.2973706) q[0];
sx q[0];
rz(-0.64592823) q[0];
sx q[0];
rz(2.5883664) q[0];
rz(-0.79752254) q[1];
sx q[1];
rz(-2.1452466) q[1];
sx q[1];
rz(-1.9901989) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0259375) q[0];
sx q[0];
rz(-2.135072) q[0];
sx q[0];
rz(-1.7274117) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7417479) q[2];
sx q[2];
rz(-2.1151849) q[2];
sx q[2];
rz(0.65786568) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0935191) q[1];
sx q[1];
rz(-2.8052605) q[1];
sx q[1];
rz(0.63365714) q[1];
rz(-pi) q[2];
rz(-1.0955515) q[3];
sx q[3];
rz(-2.4184347) q[3];
sx q[3];
rz(-2.1840546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.79661757) q[2];
sx q[2];
rz(-1.4130219) q[2];
sx q[2];
rz(2.3100992) q[2];
rz(-1.3384532) q[3];
sx q[3];
rz(-2.0667388) q[3];
sx q[3];
rz(-2.2510546) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3727342) q[0];
sx q[0];
rz(-1.2059728) q[0];
sx q[0];
rz(0.15705577) q[0];
rz(-1.318469) q[1];
sx q[1];
rz(-2.1993957) q[1];
sx q[1];
rz(-0.35776055) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45159949) q[0];
sx q[0];
rz(-2.4781961) q[0];
sx q[0];
rz(2.3729352) q[0];
rz(-pi) q[1];
rz(-0.4464726) q[2];
sx q[2];
rz(-2.1887767) q[2];
sx q[2];
rz(-0.22836049) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4688229) q[1];
sx q[1];
rz(-2.2190418) q[1];
sx q[1];
rz(-2.286377) q[1];
rz(-pi) q[2];
rz(2.5671183) q[3];
sx q[3];
rz(-1.4374497) q[3];
sx q[3];
rz(-2.3337618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99178189) q[2];
sx q[2];
rz(-1.2573743) q[2];
sx q[2];
rz(3.120976) q[2];
rz(-1.2425544) q[3];
sx q[3];
rz(-0.81762448) q[3];
sx q[3];
rz(-2.8500565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6108516) q[0];
sx q[0];
rz(-1.956097) q[0];
sx q[0];
rz(2.8148742) q[0];
rz(-0.84699455) q[1];
sx q[1];
rz(-2.0920483) q[1];
sx q[1];
rz(1.0583896) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3922472) q[0];
sx q[0];
rz(-1.5010745) q[0];
sx q[0];
rz(1.0564984) q[0];
x q[1];
rz(1.7472668) q[2];
sx q[2];
rz(-0.93910774) q[2];
sx q[2];
rz(1.2588009) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8135735) q[1];
sx q[1];
rz(-1.18888) q[1];
sx q[1];
rz(3.0573634) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31553573) q[3];
sx q[3];
rz(-1.0275176) q[3];
sx q[3];
rz(1.735294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1191001) q[2];
sx q[2];
rz(-2.7969226) q[2];
sx q[2];
rz(1.1723088) q[2];
rz(-0.0017702866) q[3];
sx q[3];
rz(-2.6383196) q[3];
sx q[3];
rz(0.9790023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86972791) q[0];
sx q[0];
rz(-0.80369049) q[0];
sx q[0];
rz(1.3386238) q[0];
rz(1.9610693) q[1];
sx q[1];
rz(-1.0931284) q[1];
sx q[1];
rz(0.48042935) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.521489) q[0];
sx q[0];
rz(-1.1166945) q[0];
sx q[0];
rz(0.2257077) q[0];
rz(1.5283723) q[2];
sx q[2];
rz(-1.7406751) q[2];
sx q[2];
rz(0.9695878) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5753382) q[1];
sx q[1];
rz(-1.6702139) q[1];
sx q[1];
rz(1.0097617) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35828405) q[3];
sx q[3];
rz(-1.3582152) q[3];
sx q[3];
rz(-0.5396809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43634811) q[2];
sx q[2];
rz(-0.99506012) q[2];
sx q[2];
rz(-2.056541) q[2];
rz(2.0294225) q[3];
sx q[3];
rz(-0.89434904) q[3];
sx q[3];
rz(-0.81356847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
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
rz(0.28806624) q[0];
sx q[0];
rz(-0.94678322) q[0];
sx q[0];
rz(2.1296401) q[0];
rz(-2.2400253) q[1];
sx q[1];
rz(-2.8755867) q[1];
sx q[1];
rz(-2.576135) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8249567) q[0];
sx q[0];
rz(-0.68616435) q[0];
sx q[0];
rz(1.9030722) q[0];
x q[1];
rz(-0.26763518) q[2];
sx q[2];
rz(-1.9861172) q[2];
sx q[2];
rz(-0.079697996) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9931075) q[1];
sx q[1];
rz(-2.4186181) q[1];
sx q[1];
rz(-1.430368) q[1];
rz(-pi) q[2];
rz(-1.8875445) q[3];
sx q[3];
rz(-1.4707139) q[3];
sx q[3];
rz(0.36857482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7908343) q[2];
sx q[2];
rz(-1.5956722) q[2];
sx q[2];
rz(1.8809543) q[2];
rz(0.44220051) q[3];
sx q[3];
rz(-2.111777) q[3];
sx q[3];
rz(1.376576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.5762536) q[0];
sx q[0];
rz(-0.88247846) q[0];
sx q[0];
rz(-0.89378617) q[0];
rz(1.5671989) q[1];
sx q[1];
rz(-1.6849453) q[1];
sx q[1];
rz(-1.4243855) q[1];
rz(0.5657351) q[2];
sx q[2];
rz(-2.6489352) q[2];
sx q[2];
rz(0.87633662) q[2];
rz(2.9775053) q[3];
sx q[3];
rz(-1.224095) q[3];
sx q[3];
rz(-1.7188354) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
