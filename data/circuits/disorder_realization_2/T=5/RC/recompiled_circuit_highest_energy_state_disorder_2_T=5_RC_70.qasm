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
rz(-1.4827363) q[0];
sx q[0];
rz(4.1229376) q[0];
sx q[0];
rz(11.468588) q[0];
rz(2.3535347) q[1];
sx q[1];
rz(5.2323879) q[1];
sx q[1];
rz(9.431737) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1084676) q[0];
sx q[0];
rz(-1.5146717) q[0];
sx q[0];
rz(-2.2214487) q[0];
rz(-pi) q[1];
rz(2.2175781) q[2];
sx q[2];
rz(-1.7597464) q[2];
sx q[2];
rz(1.3786157) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6485868) q[1];
sx q[1];
rz(-1.3351591) q[1];
sx q[1];
rz(0.034502397) q[1];
x q[2];
rz(1.4640635) q[3];
sx q[3];
rz(-1.9088703) q[3];
sx q[3];
rz(2.6298912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9030582) q[2];
sx q[2];
rz(-1.7944585) q[2];
sx q[2];
rz(1.6713589) q[2];
rz(-3.0155731) q[3];
sx q[3];
rz(-2.6991548) q[3];
sx q[3];
rz(0.68388763) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11494342) q[0];
sx q[0];
rz(-0.4489972) q[0];
sx q[0];
rz(0.63585109) q[0];
rz(-0.77330971) q[1];
sx q[1];
rz(-2.4644303) q[1];
sx q[1];
rz(1.4453452) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6757386) q[0];
sx q[0];
rz(-2.2096118) q[0];
sx q[0];
rz(2.5771228) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3582621) q[2];
sx q[2];
rz(-1.223067) q[2];
sx q[2];
rz(-1.7933947) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8938039) q[1];
sx q[1];
rz(-0.87444982) q[1];
sx q[1];
rz(-1.7764093) q[1];
rz(3.0731455) q[3];
sx q[3];
rz(-1.7020149) q[3];
sx q[3];
rz(2.69773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.14629743) q[2];
sx q[2];
rz(-1.3714906) q[2];
sx q[2];
rz(-3.0038317) q[2];
rz(0.45608258) q[3];
sx q[3];
rz(-2.5465953) q[3];
sx q[3];
rz(-0.47541398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59738961) q[0];
sx q[0];
rz(-1.5837639) q[0];
sx q[0];
rz(0.57918817) q[0];
rz(3.0384565) q[1];
sx q[1];
rz(-1.8201647) q[1];
sx q[1];
rz(0.78027049) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53952459) q[0];
sx q[0];
rz(-1.4756265) q[0];
sx q[0];
rz(1.5896941) q[0];
rz(-pi) q[1];
rz(-1.7306855) q[2];
sx q[2];
rz(-0.89582755) q[2];
sx q[2];
rz(2.8596148) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0611813) q[1];
sx q[1];
rz(-1.1763417) q[1];
sx q[1];
rz(0.2385255) q[1];
rz(-2.2663651) q[3];
sx q[3];
rz(-1.4000633) q[3];
sx q[3];
rz(0.010893498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66037336) q[2];
sx q[2];
rz(-1.4619091) q[2];
sx q[2];
rz(0.090864651) q[2];
rz(-1.6615435) q[3];
sx q[3];
rz(-1.3250947) q[3];
sx q[3];
rz(0.81502325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12260967) q[0];
sx q[0];
rz(-2.9538587) q[0];
sx q[0];
rz(2.6222099) q[0];
rz(1.9620365) q[1];
sx q[1];
rz(-2.4603381) q[1];
sx q[1];
rz(0.048351668) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78084842) q[0];
sx q[0];
rz(-2.0791302) q[0];
sx q[0];
rz(-0.8698277) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7555356) q[2];
sx q[2];
rz(-2.0252844) q[2];
sx q[2];
rz(-2.7234361) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0569339) q[1];
sx q[1];
rz(-2.6531583) q[1];
sx q[1];
rz(-3.0822251) q[1];
x q[2];
rz(-2.3465041) q[3];
sx q[3];
rz(-1.7582736) q[3];
sx q[3];
rz(-0.80611967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4590596) q[2];
sx q[2];
rz(-2.8352663) q[2];
sx q[2];
rz(1.4502067) q[2];
rz(2.0356483) q[3];
sx q[3];
rz(-1.148843) q[3];
sx q[3];
rz(2.2753184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-0.76097101) q[0];
sx q[0];
rz(-2.3411317) q[0];
sx q[0];
rz(2.3642484) q[0];
rz(-2.6457973) q[1];
sx q[1];
rz(-0.40758857) q[1];
sx q[1];
rz(0.19283238) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58932006) q[0];
sx q[0];
rz(-2.5077132) q[0];
sx q[0];
rz(0.53323563) q[0];
rz(-pi) q[1];
rz(-0.2482581) q[2];
sx q[2];
rz(-0.78310668) q[2];
sx q[2];
rz(-0.017489028) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7954171) q[1];
sx q[1];
rz(-0.3095135) q[1];
sx q[1];
rz(-1.8619821) q[1];
rz(-pi) q[2];
rz(0.16280547) q[3];
sx q[3];
rz(-1.9678402) q[3];
sx q[3];
rz(-2.5474078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5567646) q[2];
sx q[2];
rz(-0.77797055) q[2];
sx q[2];
rz(-0.83621109) q[2];
rz(0.95544514) q[3];
sx q[3];
rz(-1.7183869) q[3];
sx q[3];
rz(-0.53171617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8504976) q[0];
sx q[0];
rz(-1.2245155) q[0];
sx q[0];
rz(3.1078872) q[0];
rz(-2.2233502) q[1];
sx q[1];
rz(-0.70437175) q[1];
sx q[1];
rz(-2.4423626) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0110737) q[0];
sx q[0];
rz(-1.0193745) q[0];
sx q[0];
rz(-2.5370777) q[0];
rz(-pi) q[1];
rz(0.18728687) q[2];
sx q[2];
rz(-1.5178871) q[2];
sx q[2];
rz(2.1250181) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.81046178) q[1];
sx q[1];
rz(-1.7469566) q[1];
sx q[1];
rz(2.6947652) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27856234) q[3];
sx q[3];
rz(-1.8980366) q[3];
sx q[3];
rz(0.95765169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0596727) q[2];
sx q[2];
rz(-2.4884188) q[2];
sx q[2];
rz(-0.095452249) q[2];
rz(-1.0517612) q[3];
sx q[3];
rz(-1.2442518) q[3];
sx q[3];
rz(-2.2473647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28602257) q[0];
sx q[0];
rz(-2.6595071) q[0];
sx q[0];
rz(-1.1676189) q[0];
rz(0.35320148) q[1];
sx q[1];
rz(-2.628852) q[1];
sx q[1];
rz(2.8599427) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62171157) q[0];
sx q[0];
rz(-0.48292749) q[0];
sx q[0];
rz(1.3819225) q[0];
rz(-pi) q[1];
rz(-1.4640305) q[2];
sx q[2];
rz(-2.2153478) q[2];
sx q[2];
rz(1.0463971) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0822009) q[1];
sx q[1];
rz(-1.4363465) q[1];
sx q[1];
rz(-1.7198635) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1280888) q[3];
sx q[3];
rz(-1.1508599) q[3];
sx q[3];
rz(2.1518663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1436651) q[2];
sx q[2];
rz(-1.2219656) q[2];
sx q[2];
rz(-1.3227051) q[2];
rz(-1.5023242) q[3];
sx q[3];
rz(-1.084525) q[3];
sx q[3];
rz(0.098202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99763501) q[0];
sx q[0];
rz(-1.8956381) q[0];
sx q[0];
rz(-0.27398807) q[0];
rz(-0.59448376) q[1];
sx q[1];
rz(-1.7285873) q[1];
sx q[1];
rz(0.53057539) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45121058) q[0];
sx q[0];
rz(-1.57461) q[0];
sx q[0];
rz(-0.03937748) q[0];
rz(-pi) q[1];
rz(-0.55270393) q[2];
sx q[2];
rz(-1.505902) q[2];
sx q[2];
rz(-0.12285025) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.12996) q[1];
sx q[1];
rz(-1.7372565) q[1];
sx q[1];
rz(-0.37977207) q[1];
rz(1.7591997) q[3];
sx q[3];
rz(-0.72160463) q[3];
sx q[3];
rz(2.5327275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8990367) q[2];
sx q[2];
rz(-1.9344923) q[2];
sx q[2];
rz(0.25699082) q[2];
rz(1.1993923) q[3];
sx q[3];
rz(-1.2446087) q[3];
sx q[3];
rz(-1.5132343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5892107) q[0];
sx q[0];
rz(-2.1871545) q[0];
sx q[0];
rz(-2.6115665) q[0];
rz(-1.5026622) q[1];
sx q[1];
rz(-1.9866147) q[1];
sx q[1];
rz(0.93856215) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6185222) q[0];
sx q[0];
rz(-1.1022727) q[0];
sx q[0];
rz(0.65798385) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.45088) q[2];
sx q[2];
rz(-1.4282104) q[2];
sx q[2];
rz(0.8220807) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5512498) q[1];
sx q[1];
rz(-1.9971202) q[1];
sx q[1];
rz(-0.23386441) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13646941) q[3];
sx q[3];
rz(-2.1519063) q[3];
sx q[3];
rz(2.7032397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.75371257) q[2];
sx q[2];
rz(-0.91999274) q[2];
sx q[2];
rz(2.06125) q[2];
rz(2.3267817) q[3];
sx q[3];
rz(-1.1956513) q[3];
sx q[3];
rz(-1.7241534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1431047) q[0];
sx q[0];
rz(-1.4838706) q[0];
sx q[0];
rz(-1.0888354) q[0];
rz(-3.0740652) q[1];
sx q[1];
rz(-1.8636401) q[1];
sx q[1];
rz(0.75928226) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7831206) q[0];
sx q[0];
rz(-1.4519339) q[0];
sx q[0];
rz(1.8274183) q[0];
rz(-pi) q[1];
rz(-1.6948872) q[2];
sx q[2];
rz(-1.3465836) q[2];
sx q[2];
rz(-1.8607163) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6435561) q[1];
sx q[1];
rz(-0.78460556) q[1];
sx q[1];
rz(-2.8348009) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57709007) q[3];
sx q[3];
rz(-2.7365757) q[3];
sx q[3];
rz(-1.4228203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33409432) q[2];
sx q[2];
rz(-1.0958025) q[2];
sx q[2];
rz(-2.5661772) q[2];
rz(0.56257644) q[3];
sx q[3];
rz(-2.176599) q[3];
sx q[3];
rz(1.7687198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060487735) q[0];
sx q[0];
rz(-1.8555547) q[0];
sx q[0];
rz(-0.9077358) q[0];
rz(-1.0870712) q[1];
sx q[1];
rz(-2.0094951) q[1];
sx q[1];
rz(3.1241945) q[1];
rz(2.9870177) q[2];
sx q[2];
rz(-1.4308725) q[2];
sx q[2];
rz(1.9785656) q[2];
rz(-0.51856507) q[3];
sx q[3];
rz(-0.66558481) q[3];
sx q[3];
rz(-0.44105327) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
