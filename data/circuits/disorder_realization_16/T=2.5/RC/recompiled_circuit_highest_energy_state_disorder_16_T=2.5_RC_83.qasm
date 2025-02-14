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
rz(0.056869153) q[0];
sx q[0];
rz(2.9480204) q[0];
sx q[0];
rz(9.8326346) q[0];
rz(-0.017539311) q[1];
sx q[1];
rz(-1.8899625) q[1];
sx q[1];
rz(-1.6012021) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6787036) q[0];
sx q[0];
rz(-2.7015831) q[0];
sx q[0];
rz(-1.4567039) q[0];
rz(-pi) q[1];
rz(0.27215927) q[2];
sx q[2];
rz(-0.46900422) q[2];
sx q[2];
rz(0.16527612) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9484798) q[1];
sx q[1];
rz(-0.44572898) q[1];
sx q[1];
rz(-0.50909252) q[1];
rz(0.82117053) q[3];
sx q[3];
rz(-1.0945012) q[3];
sx q[3];
rz(-1.1007635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5539598) q[2];
sx q[2];
rz(-1.0593869) q[2];
sx q[2];
rz(2.9239192) q[2];
rz(-2.830128) q[3];
sx q[3];
rz(-0.61190999) q[3];
sx q[3];
rz(1.9757087) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72164732) q[0];
sx q[0];
rz(-2.8657275) q[0];
sx q[0];
rz(-2.4496147) q[0];
rz(-2.3852589) q[1];
sx q[1];
rz(-0.5223918) q[1];
sx q[1];
rz(-1.5914894) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67705757) q[0];
sx q[0];
rz(-1.4486396) q[0];
sx q[0];
rz(2.4273901) q[0];
rz(-pi) q[1];
rz(1.1728376) q[2];
sx q[2];
rz(-1.474547) q[2];
sx q[2];
rz(2.725696) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8263232) q[1];
sx q[1];
rz(-0.20139748) q[1];
sx q[1];
rz(-0.65314318) q[1];
rz(0.85568537) q[3];
sx q[3];
rz(-0.69310235) q[3];
sx q[3];
rz(-2.2732796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1456566) q[2];
sx q[2];
rz(-0.17309509) q[2];
sx q[2];
rz(0.23925979) q[2];
rz(-1.4907106) q[3];
sx q[3];
rz(-1.8473293) q[3];
sx q[3];
rz(1.9132805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2759129) q[0];
sx q[0];
rz(-0.90621197) q[0];
sx q[0];
rz(-3.1269585) q[0];
rz(0.89302653) q[1];
sx q[1];
rz(-2.1664186) q[1];
sx q[1];
rz(-0.21569529) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1050284) q[0];
sx q[0];
rz(-1.6220548) q[0];
sx q[0];
rz(3.1369484) q[0];
rz(-1.2273618) q[2];
sx q[2];
rz(-2.2015155) q[2];
sx q[2];
rz(1.3892) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0552935) q[1];
sx q[1];
rz(-2.8759416) q[1];
sx q[1];
rz(-2.6363027) q[1];
x q[2];
rz(2.3560432) q[3];
sx q[3];
rz(-2.2758898) q[3];
sx q[3];
rz(2.7712703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0821685) q[2];
sx q[2];
rz(-1.9614204) q[2];
sx q[2];
rz(-2.1480985) q[2];
rz(-2.4332186) q[3];
sx q[3];
rz(-2.0230484) q[3];
sx q[3];
rz(-0.12349252) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3941536) q[0];
sx q[0];
rz(-2.6469632) q[0];
sx q[0];
rz(-2.4476449) q[0];
rz(-2.8548062) q[1];
sx q[1];
rz(-1.6096121) q[1];
sx q[1];
rz(-2.0538816) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74831731) q[0];
sx q[0];
rz(-1.4469742) q[0];
sx q[0];
rz(-1.0626777) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7544341) q[2];
sx q[2];
rz(-0.70068073) q[2];
sx q[2];
rz(0.59631077) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3301311) q[1];
sx q[1];
rz(-2.4566659) q[1];
sx q[1];
rz(1.0792988) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4272653) q[3];
sx q[3];
rz(-1.5363524) q[3];
sx q[3];
rz(1.2726651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9168758) q[2];
sx q[2];
rz(-1.1563053) q[2];
sx q[2];
rz(-1.7737596) q[2];
rz(-1.2380838) q[3];
sx q[3];
rz(-1.5725458) q[3];
sx q[3];
rz(2.9253173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.5991768) q[0];
sx q[0];
rz(-1.8734064) q[0];
sx q[0];
rz(-0.61781484) q[0];
rz(-1.1082331) q[1];
sx q[1];
rz(-2.8384659) q[1];
sx q[1];
rz(-0.88315001) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9734479) q[0];
sx q[0];
rz(-0.75605481) q[0];
sx q[0];
rz(3.130828) q[0];
rz(-2.8034535) q[2];
sx q[2];
rz(-2.6827742) q[2];
sx q[2];
rz(0.79566075) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0269491) q[1];
sx q[1];
rz(-2.0609239) q[1];
sx q[1];
rz(3.0672706) q[1];
rz(0.41188052) q[3];
sx q[3];
rz(-1.05384) q[3];
sx q[3];
rz(-0.98983562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9250179) q[2];
sx q[2];
rz(-2.8897132) q[2];
sx q[2];
rz(1.3612755) q[2];
rz(-1.8682293) q[3];
sx q[3];
rz(-1.4342156) q[3];
sx q[3];
rz(-2.1586965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0971766) q[0];
sx q[0];
rz(-1.8955078) q[0];
sx q[0];
rz(-2.521305) q[0];
rz(-0.39707956) q[1];
sx q[1];
rz(-0.47080165) q[1];
sx q[1];
rz(-1.1995859) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7714149) q[0];
sx q[0];
rz(-1.6770898) q[0];
sx q[0];
rz(3.1205157) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.994147) q[2];
sx q[2];
rz(-1.7803528) q[2];
sx q[2];
rz(-0.042590754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7677557) q[1];
sx q[1];
rz(-1.3779614) q[1];
sx q[1];
rz(2.6111433) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34440094) q[3];
sx q[3];
rz(-1.6823497) q[3];
sx q[3];
rz(-1.2408181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24641307) q[2];
sx q[2];
rz(-1.2440888) q[2];
sx q[2];
rz(-0.67071521) q[2];
rz(2.8504168) q[3];
sx q[3];
rz(-0.61101919) q[3];
sx q[3];
rz(-2.5197855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7156242) q[0];
sx q[0];
rz(-1.427303) q[0];
sx q[0];
rz(-2.9631462) q[0];
rz(0.87002358) q[1];
sx q[1];
rz(-2.2585637) q[1];
sx q[1];
rz(2.3387486) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3960796) q[0];
sx q[0];
rz(-0.95142309) q[0];
sx q[0];
rz(-1.6016763) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8151547) q[2];
sx q[2];
rz(-2.1409224) q[2];
sx q[2];
rz(1.9889835) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.57299388) q[1];
sx q[1];
rz(-0.90403176) q[1];
sx q[1];
rz(1.9896889) q[1];
rz(-pi) q[2];
rz(2.1986897) q[3];
sx q[3];
rz(-1.3602837) q[3];
sx q[3];
rz(-0.26417662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.37658438) q[2];
sx q[2];
rz(-1.2976126) q[2];
sx q[2];
rz(0.89361781) q[2];
rz(-1.5997959) q[3];
sx q[3];
rz(-1.6518281) q[3];
sx q[3];
rz(-2.5347575) q[3];
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
rz(2.1641418) q[0];
sx q[0];
rz(-2.7315388) q[0];
sx q[0];
rz(-0.2151016) q[0];
rz(0.84456259) q[1];
sx q[1];
rz(-1.3203878) q[1];
sx q[1];
rz(0.79427687) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71090224) q[0];
sx q[0];
rz(-1.5973418) q[0];
sx q[0];
rz(1.5826167) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.568001) q[2];
sx q[2];
rz(-2.6124138) q[2];
sx q[2];
rz(-2.0647788) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6805598) q[1];
sx q[1];
rz(-0.54745142) q[1];
sx q[1];
rz(-1.3957109) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32799649) q[3];
sx q[3];
rz(-0.80910002) q[3];
sx q[3];
rz(-2.3878333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.27791417) q[2];
sx q[2];
rz(-1.2028376) q[2];
sx q[2];
rz(1.7318783) q[2];
rz(-2.388741) q[3];
sx q[3];
rz(-1.1466305) q[3];
sx q[3];
rz(-2.1394155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79733217) q[0];
sx q[0];
rz(-2.7334038) q[0];
sx q[0];
rz(2.1687188) q[0];
rz(-1.5219888) q[1];
sx q[1];
rz(-1.3643967) q[1];
sx q[1];
rz(-0.63046986) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75062932) q[0];
sx q[0];
rz(-1.4328151) q[0];
sx q[0];
rz(2.9580381) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0075118709) q[2];
sx q[2];
rz(-1.7277311) q[2];
sx q[2];
rz(-0.73534009) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9520397) q[1];
sx q[1];
rz(-1.6206883) q[1];
sx q[1];
rz(0.76811172) q[1];
rz(-2.7961066) q[3];
sx q[3];
rz(-0.18426963) q[3];
sx q[3];
rz(-0.94217448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1260599) q[2];
sx q[2];
rz(-1.9893179) q[2];
sx q[2];
rz(1.3014334) q[2];
rz(1.5850916) q[3];
sx q[3];
rz(-2.3157178) q[3];
sx q[3];
rz(2.7263156) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5179317) q[0];
sx q[0];
rz(-2.9908337) q[0];
sx q[0];
rz(2.6483722) q[0];
rz(1.2552235) q[1];
sx q[1];
rz(-1.8938277) q[1];
sx q[1];
rz(-0.36466041) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4987558) q[0];
sx q[0];
rz(-2.7793573) q[0];
sx q[0];
rz(1.1722159) q[0];
rz(-pi) q[1];
rz(-0.45367806) q[2];
sx q[2];
rz(-2.1725906) q[2];
sx q[2];
rz(-1.2860677) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9182049) q[1];
sx q[1];
rz(-2.3406696) q[1];
sx q[1];
rz(-0.078322874) q[1];
rz(-3.10121) q[3];
sx q[3];
rz(-1.2546439) q[3];
sx q[3];
rz(-2.169956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4517335) q[2];
sx q[2];
rz(-1.006459) q[2];
sx q[2];
rz(-2.149392) q[2];
rz(-2.774488) q[3];
sx q[3];
rz(-1.4398984) q[3];
sx q[3];
rz(-0.67830694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.622396) q[0];
sx q[0];
rz(-1.4764897) q[0];
sx q[0];
rz(1.3523703) q[0];
rz(-2.2182111) q[1];
sx q[1];
rz(-0.8538178) q[1];
sx q[1];
rz(-1.1710844) q[1];
rz(0.14262234) q[2];
sx q[2];
rz(-1.3726948) q[2];
sx q[2];
rz(1.2683327) q[2];
rz(0.67305858) q[3];
sx q[3];
rz(-2.06235) q[3];
sx q[3];
rz(1.6497117) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
