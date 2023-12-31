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
rz(-1.4562343) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0460912) q[0];
sx q[0];
rz(-1.7016194) q[0];
sx q[0];
rz(-0.019898947) q[0];
x q[1];
rz(0.37402447) q[2];
sx q[2];
rz(-2.0846539) q[2];
sx q[2];
rz(0.49486578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.734182) q[1];
sx q[1];
rz(-1.721518) q[1];
sx q[1];
rz(0.68024866) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9772894) q[3];
sx q[3];
rz(-0.32391732) q[3];
sx q[3];
rz(1.2795554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.75498092) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(2.6888729) q[2];
rz(0.1581986) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(-2.246777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2125856) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(-1.1520977) q[0];
rz(1.903803) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(0.47098413) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1176227) q[0];
sx q[0];
rz(-1.6256486) q[0];
sx q[0];
rz(-2.0532002) q[0];
x q[1];
rz(-0.97371308) q[2];
sx q[2];
rz(-1.7657585) q[2];
sx q[2];
rz(0.82414579) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.7324595) q[1];
sx q[1];
rz(-2.0583378) q[1];
sx q[1];
rz(0.061846102) q[1];
rz(2.4460375) q[3];
sx q[3];
rz(-1.5091981) q[3];
sx q[3];
rz(1.0775281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0835138) q[2];
sx q[2];
rz(-0.60516548) q[2];
sx q[2];
rz(2.8857968) q[2];
rz(-1.4852218) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(-2.9799057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8783022) q[0];
sx q[0];
rz(-2.0354164) q[0];
sx q[0];
rz(0.27134744) q[0];
rz(-0.73633206) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(2.7022865) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96734756) q[0];
sx q[0];
rz(-0.085701533) q[0];
sx q[0];
rz(1.218319) q[0];
x q[1];
rz(-1.5191684) q[2];
sx q[2];
rz(-1.2106967) q[2];
sx q[2];
rz(2.1215631) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9648793) q[1];
sx q[1];
rz(-0.63767725) q[1];
sx q[1];
rz(-1.2769075) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2112591) q[3];
sx q[3];
rz(-0.6593245) q[3];
sx q[3];
rz(-1.7742771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.97418857) q[2];
sx q[2];
rz(-1.5891275) q[2];
sx q[2];
rz(2.9411194) q[2];
rz(2.386507) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(0.42373207) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0641091) q[0];
sx q[0];
rz(-1.5492726) q[0];
sx q[0];
rz(1.8970998) q[0];
rz(-2.3311133) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(2.2669852) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4865206) q[0];
sx q[0];
rz(-2.1694896) q[0];
sx q[0];
rz(-0.18874164) q[0];
x q[1];
rz(0.34747296) q[2];
sx q[2];
rz(-2.7345737) q[2];
sx q[2];
rz(2.4328872) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7728459) q[1];
sx q[1];
rz(-2.6987942) q[1];
sx q[1];
rz(-0.55261353) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0365385) q[3];
sx q[3];
rz(-1.6451562) q[3];
sx q[3];
rz(-0.33945938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.13742927) q[2];
sx q[2];
rz(-0.88976294) q[2];
sx q[2];
rz(-1.4271663) q[2];
rz(0.066120474) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(-1.6920413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78113294) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(-2.5710035) q[0];
rz(0.55496201) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(-0.74329174) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3742204) q[0];
sx q[0];
rz(-1.351007) q[0];
sx q[0];
rz(-0.90669294) q[0];
x q[1];
rz(1.6749803) q[2];
sx q[2];
rz(-1.0717234) q[2];
sx q[2];
rz(-0.68945976) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44504657) q[1];
sx q[1];
rz(-1.7940709) q[1];
sx q[1];
rz(2.9162507) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27637847) q[3];
sx q[3];
rz(-0.82377269) q[3];
sx q[3];
rz(0.36229047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9479998) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(-2.8732079) q[2];
rz(2.0949481) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(2.823765) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.627581) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(-0.50672379) q[0];
rz(-0.2535893) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(-2.0862897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9133638) q[0];
sx q[0];
rz(-2.1852487) q[0];
sx q[0];
rz(0.22426228) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52168092) q[2];
sx q[2];
rz(-0.98419596) q[2];
sx q[2];
rz(0.66643836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9195337) q[1];
sx q[1];
rz(-2.2026081) q[1];
sx q[1];
rz(2.6421089) q[1];
x q[2];
rz(-1.8063227) q[3];
sx q[3];
rz(-0.84170656) q[3];
sx q[3];
rz(1.0585166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48173299) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(-2.2591023) q[2];
rz(-0.64583889) q[3];
sx q[3];
rz(-1.9927988) q[3];
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
rz(-pi) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5053453) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(0.77254599) q[0];
rz(-1.729471) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(-0.57377446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6928497) q[0];
sx q[0];
rz(-0.22073711) q[0];
sx q[0];
rz(1.0462532) q[0];
rz(-1.3705809) q[2];
sx q[2];
rz(-1.5491312) q[2];
sx q[2];
rz(1.3862762) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.335865) q[1];
sx q[1];
rz(-2.1153567) q[1];
sx q[1];
rz(-0.90421275) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9912234) q[3];
sx q[3];
rz(-1.8024076) q[3];
sx q[3];
rz(1.3814955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1022169) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(0.77073628) q[2];
rz(-0.43631521) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(1.2873945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
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
rz(-1.8687826) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(1.9708721) q[0];
rz(-2.6314578) q[1];
sx q[1];
rz(-1.7850103) q[1];
sx q[1];
rz(1.8458813) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1575748) q[0];
sx q[0];
rz(-0.57849738) q[0];
sx q[0];
rz(2.5111141) q[0];
x q[1];
rz(-2.2175118) q[2];
sx q[2];
rz(-2.7087822) q[2];
sx q[2];
rz(-0.88976394) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8746652) q[1];
sx q[1];
rz(-1.7942567) q[1];
sx q[1];
rz(0.93519559) q[1];
rz(-pi) q[2];
rz(-2.0115764) q[3];
sx q[3];
rz(-1.545558) q[3];
sx q[3];
rz(2.4079635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7065113) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(-0.56813017) q[2];
rz(2.1614697) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(-2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(-1.0124857) q[1];
sx q[1];
rz(0.83818865) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4351589) q[0];
sx q[0];
rz(-2.7088232) q[0];
sx q[0];
rz(-0.5612527) q[0];
rz(-pi) q[1];
rz(-0.067473472) q[2];
sx q[2];
rz(-0.81574342) q[2];
sx q[2];
rz(0.93630723) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5327685) q[1];
sx q[1];
rz(-2.3012487) q[1];
sx q[1];
rz(-2.8735012) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8348654) q[3];
sx q[3];
rz(-2.8711651) q[3];
sx q[3];
rz(-0.21817792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80205408) q[2];
sx q[2];
rz(-2.4513117) q[2];
sx q[2];
rz(1.8748803) q[2];
rz(0.96261111) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4203913) q[0];
sx q[0];
rz(-1.2370011) q[0];
sx q[0];
rz(-0.23751968) q[0];
rz(1.0182084) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(-0.231803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12122614) q[0];
sx q[0];
rz(-2.0775284) q[0];
sx q[0];
rz(3.0324742) q[0];
rz(-2.899029) q[2];
sx q[2];
rz(-1.3360268) q[2];
sx q[2];
rz(-2.501542) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4494891) q[1];
sx q[1];
rz(-0.30174258) q[1];
sx q[1];
rz(2.6471789) q[1];
x q[2];
rz(0.2139123) q[3];
sx q[3];
rz(-2.247346) q[3];
sx q[3];
rz(0.032534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1101749) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(0.38816372) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(2.3378519) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5647472) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(-0.96881962) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(1.8112524) q[2];
sx q[2];
rz(-1.8869055) q[2];
sx q[2];
rz(-1.4958924) q[2];
rz(1.5512636) q[3];
sx q[3];
rz(-1.7587147) q[3];
sx q[3];
rz(1.3295909) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
