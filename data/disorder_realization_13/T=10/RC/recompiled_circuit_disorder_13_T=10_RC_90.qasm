OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.17833248) q[0];
sx q[0];
rz(-1.4890716) q[0];
sx q[0];
rz(-0.89515495) q[0];
rz(2.826638) q[1];
sx q[1];
rz(-1.0840253) q[1];
sx q[1];
rz(1.6853583) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47269868) q[0];
sx q[0];
rz(-1.5905252) q[0];
sx q[0];
rz(1.701645) q[0];
rz(0.99631359) q[2];
sx q[2];
rz(-0.62553863) q[2];
sx q[2];
rz(-1.1686981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1617042) q[1];
sx q[1];
rz(-2.4474505) q[1];
sx q[1];
rz(-2.90467) q[1];
rz(-1.5159357) q[3];
sx q[3];
rz(-1.8901955) q[3];
sx q[3];
rz(1.1064135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3866117) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(2.6888729) q[2];
rz(0.1581986) q[3];
sx q[3];
rz(-0.69163624) q[3];
sx q[3];
rz(-0.89481568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92900705) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(1.989495) q[0];
rz(1.2377897) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(-0.47098413) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5841056) q[0];
sx q[0];
rz(-0.4852681) q[0];
sx q[0];
rz(-1.6886061) q[0];
x q[1];
rz(2.1678796) q[2];
sx q[2];
rz(-1.7657585) q[2];
sx q[2];
rz(-2.3174469) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.7324595) q[1];
sx q[1];
rz(-1.0832548) q[1];
sx q[1];
rz(-3.0797466) q[1];
x q[2];
rz(0.095951565) q[3];
sx q[3];
rz(-2.4437685) q[3];
sx q[3];
rz(2.5747091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.058078893) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(2.8857968) q[2];
rz(1.6563709) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(0.16168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8783022) q[0];
sx q[0];
rz(-2.0354164) q[0];
sx q[0];
rz(2.8702452) q[0];
rz(-0.73633206) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(2.7022865) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1742451) q[0];
sx q[0];
rz(-0.085701533) q[0];
sx q[0];
rz(-1.218319) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6224242) q[2];
sx q[2];
rz(-1.2106967) q[2];
sx q[2];
rz(2.1215631) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9860552) q[1];
sx q[1];
rz(-1.3974766) q[1];
sx q[1];
rz(-0.9539414) q[1];
x q[2];
rz(-0.4337173) q[3];
sx q[3];
rz(-1.0573514) q[3];
sx q[3];
rz(-2.123326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1674041) q[2];
sx q[2];
rz(-1.5891275) q[2];
sx q[2];
rz(0.20047323) q[2];
rz(-2.386507) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(-0.42373207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0641091) q[0];
sx q[0];
rz(-1.5492726) q[0];
sx q[0];
rz(1.8970998) q[0];
rz(0.81047932) q[1];
sx q[1];
rz(-1.8119718) q[1];
sx q[1];
rz(0.8746075) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8133102) q[0];
sx q[0];
rz(-2.5173442) q[0];
sx q[0];
rz(1.3024131) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38509102) q[2];
sx q[2];
rz(-1.706012) q[2];
sx q[2];
rz(1.9584292) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9676535) q[1];
sx q[1];
rz(-1.9441009) q[1];
sx q[1];
rz(-1.3268382) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61769684) q[3];
sx q[3];
rz(-0.12862895) q[3];
sx q[3];
rz(2.5240412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.13742927) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(-1.7144263) q[2];
rz(-0.066120474) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(1.6920413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(0.78113294) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(-2.5710035) q[0];
rz(-0.55496201) q[1];
sx q[1];
rz(-0.73736063) q[1];
sx q[1];
rz(2.3983009) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76737228) q[0];
sx q[0];
rz(-1.7905856) q[0];
sx q[0];
rz(0.90669294) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50136106) q[2];
sx q[2];
rz(-1.4793581) q[2];
sx q[2];
rz(-0.83133343) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.35768269) q[1];
sx q[1];
rz(-0.3158814) q[1];
sx q[1];
rz(-0.79343474) q[1];
rz(-1.284243) q[3];
sx q[3];
rz(-0.78714579) q[3];
sx q[3];
rz(-3.1084276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1935929) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(1.0466446) q[3];
sx q[3];
rz(-0.36473754) q[3];
sx q[3];
rz(-0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.627581) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(2.6348689) q[0];
rz(2.8880033) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(2.0862897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21181606) q[0];
sx q[0];
rz(-1.388071) q[0];
sx q[0];
rz(0.94434785) q[0];
rz(0.91674532) q[2];
sx q[2];
rz(-1.9987717) q[2];
sx q[2];
rz(2.5452754) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1736974) q[1];
sx q[1];
rz(-2.3579862) q[1];
sx q[1];
rz(-0.99131363) q[1];
rz(-pi) q[2];
rz(-2.8860693) q[3];
sx q[3];
rz(-0.75948411) q[3];
sx q[3];
rz(-1.7373191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48173299) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(2.2591023) q[2];
rz(2.4957538) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(1.283949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362474) q[0];
sx q[0];
rz(-1.212965) q[0];
sx q[0];
rz(-0.77254599) q[0];
rz(1.4121217) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(-0.57377446) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.505578) q[0];
sx q[0];
rz(-1.6806707) q[0];
sx q[0];
rz(1.3789603) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4622719) q[2];
sx q[2];
rz(-0.2013686) q[2];
sx q[2];
rz(-0.078171922) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.335865) q[1];
sx q[1];
rz(-2.1153567) q[1];
sx q[1];
rz(-2.2373799) q[1];
rz(1.0049099) q[3];
sx q[3];
rz(-2.8661869) q[3];
sx q[3];
rz(-0.79808455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(-2.3708564) q[2];
rz(-2.7052774) q[3];
sx q[3];
rz(-1.2687012) q[3];
sx q[3];
rz(1.2873945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8687826) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(-1.9708721) q[0];
rz(0.51013485) q[1];
sx q[1];
rz(-1.7850103) q[1];
sx q[1];
rz(1.8458813) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1575748) q[0];
sx q[0];
rz(-2.5630953) q[0];
sx q[0];
rz(2.5111141) q[0];
rz(-1.9240575) q[2];
sx q[2];
rz(-1.3152939) q[2];
sx q[2];
rz(3.0614292) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46591972) q[1];
sx q[1];
rz(-0.95341668) q[1];
sx q[1];
rz(-0.27523756) q[1];
x q[2];
rz(0.027904228) q[3];
sx q[3];
rz(-1.1301665) q[3];
sx q[3];
rz(0.84907109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7065113) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(-2.5734625) q[2];
rz(2.1614697) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(-2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4762964) q[0];
sx q[0];
rz(-0.81403533) q[0];
sx q[0];
rz(0.67767674) q[0];
rz(2.9455345) q[1];
sx q[1];
rz(-1.0124857) q[1];
sx q[1];
rz(-0.83818865) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8294551) q[0];
sx q[0];
rz(-1.2078309) q[0];
sx q[0];
rz(1.8118993) q[0];
rz(1.499275) q[2];
sx q[2];
rz(-2.3841249) q[2];
sx q[2];
rz(-2.1069991) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6088241) q[1];
sx q[1];
rz(-0.84034398) q[1];
sx q[1];
rz(0.26809147) q[1];
x q[2];
rz(-1.6543051) q[3];
sx q[3];
rz(-1.3132846) q[3];
sx q[3];
rz(2.6058692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.80205408) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.8748803) q[2];
rz(-2.1789815) q[3];
sx q[3];
rz(-1.5714785) q[3];
sx q[3];
rz(3.0468429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4203913) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(-2.904073) q[0];
rz(-1.0182084) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(-2.9097897) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7451413) q[0];
sx q[0];
rz(-1.475435) q[0];
sx q[0];
rz(-2.0800637) q[0];
rz(-1.3292153) q[2];
sx q[2];
rz(-1.3350147) q[2];
sx q[2];
rz(-2.1533522) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9634562) q[1];
sx q[1];
rz(-1.3061211) q[1];
sx q[1];
rz(1.7174277) q[1];
rz(0.88296367) q[3];
sx q[3];
rz(-1.737088) q[3];
sx q[3];
rz(-1.4030786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1101749) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(0.38816372) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(-2.3378519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5647472) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(2.172773) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(1.8112524) q[2];
sx q[2];
rz(-1.8869055) q[2];
sx q[2];
rz(-1.4958924) q[2];
rz(-1.5512636) q[3];
sx q[3];
rz(-1.3828779) q[3];
sx q[3];
rz(-1.8120017) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
