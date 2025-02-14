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
rz(0.83346811) q[0];
sx q[0];
rz(-0.8136189) q[0];
sx q[0];
rz(1.5746434) q[0];
rz(2.0180702) q[1];
sx q[1];
rz(3.968278) q[1];
sx q[1];
rz(8.9180515) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6637092) q[0];
sx q[0];
rz(-1.0179011) q[0];
sx q[0];
rz(1.5922484) q[0];
rz(-1.0584432) q[2];
sx q[2];
rz(-1.3500642) q[2];
sx q[2];
rz(-2.1519933) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.68046698) q[1];
sx q[1];
rz(-2.8790747) q[1];
sx q[1];
rz(-2.9295488) q[1];
x q[2];
rz(-3.1024974) q[3];
sx q[3];
rz(-2.2347898) q[3];
sx q[3];
rz(1.2879666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7099521) q[2];
sx q[2];
rz(-1.1974572) q[2];
sx q[2];
rz(0.54090995) q[2];
rz(-1.8139808) q[3];
sx q[3];
rz(-1.6874467) q[3];
sx q[3];
rz(2.4732164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.79880181) q[0];
sx q[0];
rz(-1.7690161) q[0];
sx q[0];
rz(-2.0846682) q[0];
rz(2.7150555) q[1];
sx q[1];
rz(-1.6273472) q[1];
sx q[1];
rz(3.1372517) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26165923) q[0];
sx q[0];
rz(-2.3525105) q[0];
sx q[0];
rz(-1.0145778) q[0];
rz(1.0356194) q[2];
sx q[2];
rz(-2.4362724) q[2];
sx q[2];
rz(0.83503228) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8480283) q[1];
sx q[1];
rz(-1.5604158) q[1];
sx q[1];
rz(1.6498098) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53496302) q[3];
sx q[3];
rz(-1.6239979) q[3];
sx q[3];
rz(-0.0004079013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9428375) q[2];
sx q[2];
rz(-1.402907) q[2];
sx q[2];
rz(2.7665561) q[2];
rz(-3.0936721) q[3];
sx q[3];
rz(-1.5791357) q[3];
sx q[3];
rz(-2.9722049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7584943) q[0];
sx q[0];
rz(-0.57621813) q[0];
sx q[0];
rz(-1.5118442) q[0];
rz(1.0668628) q[1];
sx q[1];
rz(-1.8424415) q[1];
sx q[1];
rz(2.3323434) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.209037) q[0];
sx q[0];
rz(-1.6106298) q[0];
sx q[0];
rz(-1.5341542) q[0];
rz(-pi) q[1];
rz(1.6406443) q[2];
sx q[2];
rz(-0.85416342) q[2];
sx q[2];
rz(-1.4586064) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7748562) q[1];
sx q[1];
rz(-1.3778169) q[1];
sx q[1];
rz(0.79905135) q[1];
x q[2];
rz(0.36181216) q[3];
sx q[3];
rz(-1.568095) q[3];
sx q[3];
rz(0.55015528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2201436) q[2];
sx q[2];
rz(-0.67110944) q[2];
sx q[2];
rz(-1.6371833) q[2];
rz(0.97357059) q[3];
sx q[3];
rz(-2.3928271) q[3];
sx q[3];
rz(-1.0934632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022738986) q[0];
sx q[0];
rz(-1.5819419) q[0];
sx q[0];
rz(0.36209011) q[0];
rz(-1.3478966) q[1];
sx q[1];
rz(-1.3840414) q[1];
sx q[1];
rz(2.0787584) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84447538) q[0];
sx q[0];
rz(-1.7430291) q[0];
sx q[0];
rz(0.14396859) q[0];
x q[1];
rz(1.8824018) q[2];
sx q[2];
rz(-2.5172293) q[2];
sx q[2];
rz(-1.7758689) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3058051) q[1];
sx q[1];
rz(-1.2286245) q[1];
sx q[1];
rz(2.8306116) q[1];
rz(-pi) q[2];
rz(1.3534668) q[3];
sx q[3];
rz(-1.8193075) q[3];
sx q[3];
rz(1.0729162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7245543) q[2];
sx q[2];
rz(-0.78196708) q[2];
sx q[2];
rz(-0.057272591) q[2];
rz(-2.1752684) q[3];
sx q[3];
rz(-2.1939907) q[3];
sx q[3];
rz(1.9362601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22369497) q[0];
sx q[0];
rz(-2.2530093) q[0];
sx q[0];
rz(2.5433871) q[0];
rz(-1.5679476) q[1];
sx q[1];
rz(-1.6697829) q[1];
sx q[1];
rz(-0.15984687) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2837958) q[0];
sx q[0];
rz(-2.3295053) q[0];
sx q[0];
rz(-1.3336729) q[0];
rz(-2.8627928) q[2];
sx q[2];
rz(-0.21950102) q[2];
sx q[2];
rz(0.9448828) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5789524) q[1];
sx q[1];
rz(-0.75030164) q[1];
sx q[1];
rz(2.5681488) q[1];
rz(-pi) q[2];
x q[2];
rz(1.615404) q[3];
sx q[3];
rz(-1.6470223) q[3];
sx q[3];
rz(1.7895229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15888958) q[2];
sx q[2];
rz(-0.44124678) q[2];
sx q[2];
rz(-0.59054792) q[2];
rz(2.1741137) q[3];
sx q[3];
rz(-1.5891821) q[3];
sx q[3];
rz(-1.9578741) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443053) q[0];
sx q[0];
rz(-2.4794281) q[0];
sx q[0];
rz(1.039132) q[0];
rz(0.93152535) q[1];
sx q[1];
rz(-0.68382278) q[1];
sx q[1];
rz(-1.1087803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174579) q[0];
sx q[0];
rz(-2.0724247) q[0];
sx q[0];
rz(2.5162016) q[0];
rz(-2.3310235) q[2];
sx q[2];
rz(-2.2165073) q[2];
sx q[2];
rz(3.1217074) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0239552) q[1];
sx q[1];
rz(-1.2148571) q[1];
sx q[1];
rz(0.61564095) q[1];
rz(1.0226457) q[3];
sx q[3];
rz(-1.3555058) q[3];
sx q[3];
rz(-1.5886653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.335677) q[2];
sx q[2];
rz(-2.6695965) q[2];
sx q[2];
rz(2.8373888) q[2];
rz(2.975587) q[3];
sx q[3];
rz(-1.7793572) q[3];
sx q[3];
rz(1.6148875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6840927) q[0];
sx q[0];
rz(-2.9170211) q[0];
sx q[0];
rz(-1.4373047) q[0];
rz(2.8999088) q[1];
sx q[1];
rz(-1.4983404) q[1];
sx q[1];
rz(-0.93745747) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8305377) q[0];
sx q[0];
rz(-1.3321345) q[0];
sx q[0];
rz(2.412459) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8514013) q[2];
sx q[2];
rz(-1.6001424) q[2];
sx q[2];
rz(1.2236694) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9328578) q[1];
sx q[1];
rz(-2.2003649) q[1];
sx q[1];
rz(-2.1634794) q[1];
rz(2.695153) q[3];
sx q[3];
rz(-2.1025799) q[3];
sx q[3];
rz(-1.4908294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9044056) q[2];
sx q[2];
rz(-2.0439456) q[2];
sx q[2];
rz(0.22906765) q[2];
rz(1.9704874) q[3];
sx q[3];
rz(-1.3874715) q[3];
sx q[3];
rz(-2.4598725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1218629) q[0];
sx q[0];
rz(-0.64058146) q[0];
sx q[0];
rz(-1.7280475) q[0];
rz(1.9238663) q[1];
sx q[1];
rz(-1.4378589) q[1];
sx q[1];
rz(2.6103643) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6474924) q[0];
sx q[0];
rz(-1.7652041) q[0];
sx q[0];
rz(-1.902778) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0036471) q[2];
sx q[2];
rz(-2.3414991) q[2];
sx q[2];
rz(-2.6587644) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.574305) q[1];
sx q[1];
rz(-1.131408) q[1];
sx q[1];
rz(2.4498825) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8645499) q[3];
sx q[3];
rz(-1.862502) q[3];
sx q[3];
rz(3.0066389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5765877) q[2];
sx q[2];
rz(-1.7030623) q[2];
sx q[2];
rz(2.3459404) q[2];
rz(0.82120454) q[3];
sx q[3];
rz(-2.932565) q[3];
sx q[3];
rz(2.8388099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8177556) q[0];
sx q[0];
rz(-2.9471286) q[0];
sx q[0];
rz(1.6545464) q[0];
rz(-1.6200292) q[1];
sx q[1];
rz(-1.227102) q[1];
sx q[1];
rz(2.0797268) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26041612) q[0];
sx q[0];
rz(-1.9435391) q[0];
sx q[0];
rz(2.580852) q[0];
rz(-pi) q[1];
rz(1.7369607) q[2];
sx q[2];
rz(-0.29815692) q[2];
sx q[2];
rz(2.1782137) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3848906) q[1];
sx q[1];
rz(-1.9918009) q[1];
sx q[1];
rz(2.9744991) q[1];
rz(2.0031018) q[3];
sx q[3];
rz(-2.6902373) q[3];
sx q[3];
rz(-0.26717523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1294641) q[2];
sx q[2];
rz(-1.1295372) q[2];
sx q[2];
rz(-0.10936347) q[2];
rz(-2.8520285) q[3];
sx q[3];
rz(-1.3907631) q[3];
sx q[3];
rz(2.4880828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7107723) q[0];
sx q[0];
rz(-1.9662974) q[0];
sx q[0];
rz(1.9868504) q[0];
rz(2.7180502) q[1];
sx q[1];
rz(-1.4332899) q[1];
sx q[1];
rz(-0.67798859) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2104312) q[0];
sx q[0];
rz(-0.95936459) q[0];
sx q[0];
rz(1.0125041) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26099796) q[2];
sx q[2];
rz(-1.7903062) q[2];
sx q[2];
rz(-2.0479808) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0748234) q[1];
sx q[1];
rz(-2.4500203) q[1];
sx q[1];
rz(-2.2645386) q[1];
x q[2];
rz(-0.5918233) q[3];
sx q[3];
rz(-1.6890967) q[3];
sx q[3];
rz(2.5448256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8573711) q[2];
sx q[2];
rz(-0.81615964) q[2];
sx q[2];
rz(2.5428499) q[2];
rz(1.6019999) q[3];
sx q[3];
rz(-1.2076104) q[3];
sx q[3];
rz(-2.1046751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67958647) q[0];
sx q[0];
rz(-0.34272598) q[0];
sx q[0];
rz(1.7525679) q[0];
rz(-0.012808375) q[1];
sx q[1];
rz(-0.3777596) q[1];
sx q[1];
rz(-1.6526745) q[1];
rz(-1.4228062) q[2];
sx q[2];
rz(-2.4952546) q[2];
sx q[2];
rz(-0.15202022) q[2];
rz(-1.5690702) q[3];
sx q[3];
rz(-2.5826517) q[3];
sx q[3];
rz(3.0351372) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
