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
rz(2.2464377) q[0];
rz(2.826638) q[1];
sx q[1];
rz(-1.0840253) q[1];
sx q[1];
rz(1.6853583) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47269868) q[0];
sx q[0];
rz(-1.5510674) q[0];
sx q[0];
rz(-1.701645) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1452791) q[2];
sx q[2];
rz(-2.516054) q[2];
sx q[2];
rz(1.9728945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.734182) q[1];
sx q[1];
rz(-1.721518) q[1];
sx q[1];
rz(-0.68024866) q[1];
rz(-pi) q[2];
rz(-1.625657) q[3];
sx q[3];
rz(-1.8901955) q[3];
sx q[3];
rz(2.0351792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.75498092) q[2];
sx q[2];
rz(-1.3957916) q[2];
sx q[2];
rz(-0.45271978) q[2];
rz(2.9833941) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(-0.89481568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92900705) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(-1.1520977) q[0];
rz(-1.2377897) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(0.47098413) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5574871) q[0];
sx q[0];
rz(-0.4852681) q[0];
sx q[0];
rz(-1.4529865) q[0];
rz(1.9085625) q[2];
sx q[2];
rz(-2.517189) q[2];
sx q[2];
rz(2.117346) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60103154) q[1];
sx q[1];
rz(-0.49113501) q[1];
sx q[1];
rz(-1.6868375) q[1];
rz(0.095951565) q[3];
sx q[3];
rz(-2.4437685) q[3];
sx q[3];
rz(-0.56688353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0835138) q[2];
sx q[2];
rz(-0.60516548) q[2];
sx q[2];
rz(-0.2557959) q[2];
rz(-1.4852218) q[3];
sx q[3];
rz(-1.9519613) q[3];
sx q[3];
rz(2.9799057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8783022) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(2.8702452) q[0];
rz(-0.73633206) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(-0.43930611) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96734756) q[0];
sx q[0];
rz(-0.085701533) q[0];
sx q[0];
rz(1.9232737) q[0];
rz(-pi) q[1];
rz(1.6224242) q[2];
sx q[2];
rz(-1.2106967) q[2];
sx q[2];
rz(-1.0200295) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1767133) q[1];
sx q[1];
rz(-2.5039154) q[1];
sx q[1];
rz(1.8646851) q[1];
rz(-pi) q[2];
rz(0.4337173) q[3];
sx q[3];
rz(-1.0573514) q[3];
sx q[3];
rz(2.123326) q[3];
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
rz(-0.75508562) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077483594) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(1.8970998) q[0];
rz(-0.81047932) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(0.8746075) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65507209) q[0];
sx q[0];
rz(-0.97210303) q[0];
sx q[0];
rz(2.952851) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38509102) q[2];
sx q[2];
rz(-1.4355806) q[2];
sx q[2];
rz(-1.9584292) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8352656) q[1];
sx q[1];
rz(-1.7976465) q[1];
sx q[1];
rz(-2.7579685) q[1];
rz(-2.5238958) q[3];
sx q[3];
rz(-0.12862895) q[3];
sx q[3];
rz(2.5240412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.13742927) q[2];
sx q[2];
rz(-0.88976294) q[2];
sx q[2];
rz(1.7144263) q[2];
rz(-0.066120474) q[3];
sx q[3];
rz(-2.7757006) q[3];
sx q[3];
rz(1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604597) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(2.5710035) q[0];
rz(0.55496201) q[1];
sx q[1];
rz(-0.73736063) q[1];
sx q[1];
rz(0.74329174) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3742204) q[0];
sx q[0];
rz(-1.351007) q[0];
sx q[0];
rz(-0.90669294) q[0];
rz(-pi) q[1];
rz(2.9530753) q[2];
sx q[2];
rz(-2.6326615) q[2];
sx q[2];
rz(0.9045507) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6965461) q[1];
sx q[1];
rz(-1.3475218) q[1];
sx q[1];
rz(-0.22534196) q[1];
rz(-pi) q[2];
rz(-1.8573496) q[3];
sx q[3];
rz(-0.78714579) q[3];
sx q[3];
rz(-0.033165008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1935929) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(-1.0466446) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(-0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.627581) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(-2.6348689) q[0];
rz(2.8880033) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(-2.0862897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2282288) q[0];
sx q[0];
rz(-2.1852487) q[0];
sx q[0];
rz(2.9173304) q[0];
rz(-pi) q[1];
rz(-2.2248473) q[2];
sx q[2];
rz(-1.142821) q[2];
sx q[2];
rz(0.59631729) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9678952) q[1];
sx q[1];
rz(-2.3579862) q[1];
sx q[1];
rz(2.150279) q[1];
rz(-pi) q[2];
x q[2];
rz(1.33527) q[3];
sx q[3];
rz(-2.2998861) q[3];
sx q[3];
rz(2.083076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.48173299) q[2];
sx q[2];
rz(-2.0969756) q[2];
sx q[2];
rz(-0.8824904) q[2];
rz(0.64583889) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(-1.8576436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5053453) q[0];
sx q[0];
rz(-1.212965) q[0];
sx q[0];
rz(0.77254599) q[0];
rz(1.729471) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(-2.5678182) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6360146) q[0];
sx q[0];
rz(-1.6806707) q[0];
sx q[0];
rz(-1.7626324) q[0];
rz(-pi) q[1];
x q[1];
rz(3.119486) q[2];
sx q[2];
rz(-1.3706285) q[2];
sx q[2];
rz(-0.18891639) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.335865) q[1];
sx q[1];
rz(-2.1153567) q[1];
sx q[1];
rz(2.2373799) q[1];
rz(-pi) q[2];
rz(1.3366367) q[3];
sx q[3];
rz(-1.7171211) q[3];
sx q[3];
rz(-2.9175266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-0.45414671) q[2];
sx q[2];
rz(-0.77073628) q[2];
rz(0.43631521) q[3];
sx q[3];
rz(-1.2687012) q[3];
sx q[3];
rz(-1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8687826) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(1.9708721) q[0];
rz(2.6314578) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(1.8458813) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9840178) q[0];
sx q[0];
rz(-2.5630953) q[0];
sx q[0];
rz(-2.5111141) q[0];
rz(-pi) q[1];
rz(0.92408085) q[2];
sx q[2];
rz(-0.43281049) q[2];
sx q[2];
rz(-2.2518287) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2669275) q[1];
sx q[1];
rz(-1.3473359) q[1];
sx q[1];
rz(2.2063971) q[1];
rz(1.1300163) q[3];
sx q[3];
rz(-1.5960346) q[3];
sx q[3];
rz(0.73362918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43508139) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(0.56813017) q[2];
rz(-2.1614697) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4762964) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(0.67767674) q[0];
rz(-2.9455345) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(2.303404) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3121376) q[0];
sx q[0];
rz(-1.2078309) q[0];
sx q[0];
rz(1.8118993) q[0];
rz(3.0741192) q[2];
sx q[2];
rz(-2.3258492) q[2];
sx q[2];
rz(-0.93630723) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14324489) q[1];
sx q[1];
rz(-1.3721826) q[1];
sx q[1];
rz(-0.82223383) q[1];
rz(-pi) q[2];
rz(-0.30672726) q[3];
sx q[3];
rz(-0.2704276) q[3];
sx q[3];
rz(-2.9234147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3395386) q[2];
sx q[2];
rz(-2.4513117) q[2];
sx q[2];
rz(-1.8748803) q[2];
rz(2.1789815) q[3];
sx q[3];
rz(-1.5714785) q[3];
sx q[3];
rz(-3.0468429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72120136) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(2.904073) q[0];
rz(-2.1233842) q[1];
sx q[1];
rz(-2.2924481) q[1];
sx q[1];
rz(-2.9097897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7451413) q[0];
sx q[0];
rz(-1.475435) q[0];
sx q[0];
rz(2.0800637) q[0];
x q[1];
rz(-2.3583057) q[2];
sx q[2];
rz(-2.8056393) q[2];
sx q[2];
rz(-2.96539) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.787549) q[1];
sx q[1];
rz(-1.429306) q[1];
sx q[1];
rz(0.267412) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2139123) q[3];
sx q[3];
rz(-0.89424664) q[3];
sx q[3];
rz(-3.1090581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1101749) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(-2.0533662) q[2];
rz(-2.7534289) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(-0.80374074) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5768455) q[0];
sx q[0];
rz(-1.7871465) q[0];
sx q[0];
rz(-0.47252895) q[0];
rz(2.172773) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(2.5122535) q[2];
sx q[2];
rz(-2.7468801) q[2];
sx q[2];
rz(2.3135452) q[2];
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