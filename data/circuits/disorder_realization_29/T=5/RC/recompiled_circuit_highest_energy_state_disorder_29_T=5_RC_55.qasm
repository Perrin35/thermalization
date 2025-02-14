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
rz(-0.40773243) q[0];
sx q[0];
rz(4.3507504) q[0];
sx q[0];
rz(11.103295) q[0];
rz(0.2585803) q[1];
sx q[1];
rz(-1.9039896) q[1];
sx q[1];
rz(3.0078476) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8834849) q[0];
sx q[0];
rz(-2.9843669) q[0];
sx q[0];
rz(1.6088465) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57853477) q[2];
sx q[2];
rz(-2.0548487) q[2];
sx q[2];
rz(0.59147385) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7647296) q[1];
sx q[1];
rz(-1.8718182) q[1];
sx q[1];
rz(0.55913289) q[1];
rz(-pi) q[2];
rz(-1.7236563) q[3];
sx q[3];
rz(-1.5706078) q[3];
sx q[3];
rz(-2.3589695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3818843) q[2];
sx q[2];
rz(-2.3298161) q[2];
sx q[2];
rz(-1.831057) q[2];
rz(-1.7740907) q[3];
sx q[3];
rz(-1.1279227) q[3];
sx q[3];
rz(3.1254356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97717184) q[0];
sx q[0];
rz(-1.1559957) q[0];
sx q[0];
rz(-2.1261862) q[0];
rz(-0.39335355) q[1];
sx q[1];
rz(-1.4721847) q[1];
sx q[1];
rz(2.4400087) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16879747) q[0];
sx q[0];
rz(-1.4864362) q[0];
sx q[0];
rz(-0.85524164) q[0];
rz(-1.8737707) q[2];
sx q[2];
rz(-1.5810484) q[2];
sx q[2];
rz(1.4958737) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9394221) q[1];
sx q[1];
rz(-1.8864453) q[1];
sx q[1];
rz(3.062518) q[1];
x q[2];
rz(-2.2982992) q[3];
sx q[3];
rz(-0.95423428) q[3];
sx q[3];
rz(0.17184251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1141438) q[2];
sx q[2];
rz(-1.0237209) q[2];
sx q[2];
rz(2.0360428) q[2];
rz(2.8607199) q[3];
sx q[3];
rz(-2.5307145) q[3];
sx q[3];
rz(-2.9893379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1059145) q[0];
sx q[0];
rz(-2.4690101) q[0];
sx q[0];
rz(1.0726844) q[0];
rz(3.1283123) q[1];
sx q[1];
rz(-1.599879) q[1];
sx q[1];
rz(0.57659155) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.018533) q[0];
sx q[0];
rz(-1.6356688) q[0];
sx q[0];
rz(0.28831248) q[0];
x q[1];
rz(2.1303612) q[2];
sx q[2];
rz(-2.9850327) q[2];
sx q[2];
rz(2.8383534) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8416764) q[1];
sx q[1];
rz(-0.51989976) q[1];
sx q[1];
rz(1.9701881) q[1];
rz(1.0131888) q[3];
sx q[3];
rz(-1.3609145) q[3];
sx q[3];
rz(-0.55279532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2836634) q[2];
sx q[2];
rz(-0.89202213) q[2];
sx q[2];
rz(-0.64209765) q[2];
rz(-0.56504956) q[3];
sx q[3];
rz(-2.8665906) q[3];
sx q[3];
rz(-2.9716085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2791486) q[0];
sx q[0];
rz(-1.886241) q[0];
sx q[0];
rz(0.24236648) q[0];
rz(-0.53388059) q[1];
sx q[1];
rz(-0.68478525) q[1];
sx q[1];
rz(-1.6530316) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72172644) q[0];
sx q[0];
rz(-0.97630608) q[0];
sx q[0];
rz(1.6564903) q[0];
x q[1];
rz(1.2514417) q[2];
sx q[2];
rz(-1.1612301) q[2];
sx q[2];
rz(-1.7920272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7515392) q[1];
sx q[1];
rz(-1.281847) q[1];
sx q[1];
rz(-0.76559905) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8418301) q[3];
sx q[3];
rz(-1.5325755) q[3];
sx q[3];
rz(2.770407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36294595) q[2];
sx q[2];
rz(-1.2641509) q[2];
sx q[2];
rz(-1.6710336) q[2];
rz(3.0342024) q[3];
sx q[3];
rz(-1.3067747) q[3];
sx q[3];
rz(1.2764527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.030468) q[0];
sx q[0];
rz(-0.46115369) q[0];
sx q[0];
rz(1.2029458) q[0];
rz(2.2466834) q[1];
sx q[1];
rz(-1.1824965) q[1];
sx q[1];
rz(2.9295909) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0184191) q[0];
sx q[0];
rz(-0.97476649) q[0];
sx q[0];
rz(2.6072748) q[0];
x q[1];
rz(2.6486829) q[2];
sx q[2];
rz(-2.5625336) q[2];
sx q[2];
rz(2.2185203) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.13276671) q[1];
sx q[1];
rz(-1.0100875) q[1];
sx q[1];
rz(-0.53489034) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5552484) q[3];
sx q[3];
rz(-0.55335303) q[3];
sx q[3];
rz(-1.4912332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2682858) q[2];
sx q[2];
rz(-2.3640859) q[2];
sx q[2];
rz(3.0262465) q[2];
rz(-2.1558732) q[3];
sx q[3];
rz(-1.7306354) q[3];
sx q[3];
rz(-2.7058097) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2301521) q[0];
sx q[0];
rz(-1.6807012) q[0];
sx q[0];
rz(0.69754115) q[0];
rz(-0.41360924) q[1];
sx q[1];
rz(-1.4613084) q[1];
sx q[1];
rz(-0.70346171) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1915775) q[0];
sx q[0];
rz(-0.30739014) q[0];
sx q[0];
rz(2.6527575) q[0];
x q[1];
rz(1.037961) q[2];
sx q[2];
rz(-2.46358) q[2];
sx q[2];
rz(0.17282669) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5825962) q[1];
sx q[1];
rz(-1.1038053) q[1];
sx q[1];
rz(1.9860616) q[1];
rz(-pi) q[2];
x q[2];
rz(2.260842) q[3];
sx q[3];
rz(-0.38057571) q[3];
sx q[3];
rz(0.56169034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0387663) q[2];
sx q[2];
rz(-0.26831728) q[2];
sx q[2];
rz(1.212567) q[2];
rz(-2.8596527) q[3];
sx q[3];
rz(-1.5423256) q[3];
sx q[3];
rz(0.50659242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7956227) q[0];
sx q[0];
rz(-0.75241929) q[0];
sx q[0];
rz(0.84683007) q[0];
rz(-0.94547358) q[1];
sx q[1];
rz(-1.5739601) q[1];
sx q[1];
rz(0.6792773) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6300129) q[0];
sx q[0];
rz(-2.0468477) q[0];
sx q[0];
rz(1.6569754) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8273583) q[2];
sx q[2];
rz(-0.96254331) q[2];
sx q[2];
rz(-1.7781374) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0019466) q[1];
sx q[1];
rz(-2.0813542) q[1];
sx q[1];
rz(-0.21378784) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6369989) q[3];
sx q[3];
rz(-0.28488628) q[3];
sx q[3];
rz(1.5004304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.88963738) q[2];
sx q[2];
rz(-2.1869982) q[2];
sx q[2];
rz(1.6471222) q[2];
rz(-2.5877118) q[3];
sx q[3];
rz(-1.5364105) q[3];
sx q[3];
rz(1.4166098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81903356) q[0];
sx q[0];
rz(-0.77545866) q[0];
sx q[0];
rz(-2.0490647) q[0];
rz(-2.20772) q[1];
sx q[1];
rz(-1.7364419) q[1];
sx q[1];
rz(0.78528231) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063443156) q[0];
sx q[0];
rz(-2.5012996) q[0];
sx q[0];
rz(-0.75115738) q[0];
x q[1];
rz(-1.0243591) q[2];
sx q[2];
rz(-1.1866202) q[2];
sx q[2];
rz(-2.667493) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8173302) q[1];
sx q[1];
rz(-1.8166537) q[1];
sx q[1];
rz(-2.861633) q[1];
x q[2];
rz(-0.58061751) q[3];
sx q[3];
rz(-2.436565) q[3];
sx q[3];
rz(0.50268307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.87172047) q[2];
sx q[2];
rz(-1.4464658) q[2];
sx q[2];
rz(-0.63228697) q[2];
rz(-2.2271683) q[3];
sx q[3];
rz(-1.2605896) q[3];
sx q[3];
rz(-2.9818592) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81894994) q[0];
sx q[0];
rz(-0.85700789) q[0];
sx q[0];
rz(2.6981165) q[0];
rz(-2.8387198) q[1];
sx q[1];
rz(-2.5587406) q[1];
sx q[1];
rz(1.3076967) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1297559) q[0];
sx q[0];
rz(-2.7934953) q[0];
sx q[0];
rz(-0.51842095) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6340387) q[2];
sx q[2];
rz(-1.2967464) q[2];
sx q[2];
rz(-0.7407032) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.16567) q[1];
sx q[1];
rz(-0.56708401) q[1];
sx q[1];
rz(-1.0600863) q[1];
rz(0.29662432) q[3];
sx q[3];
rz(-2.0231749) q[3];
sx q[3];
rz(-1.5993702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0206454) q[2];
sx q[2];
rz(-2.1006613) q[2];
sx q[2];
rz(2.5282703) q[2];
rz(0.81765085) q[3];
sx q[3];
rz(-2.652467) q[3];
sx q[3];
rz(-0.31942719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15449512) q[0];
sx q[0];
rz(-1.3767865) q[0];
sx q[0];
rz(-0.35926551) q[0];
rz(2.477395) q[1];
sx q[1];
rz(-1.3281053) q[1];
sx q[1];
rz(2.7204303) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65972787) q[0];
sx q[0];
rz(-0.43584985) q[0];
sx q[0];
rz(-1.0208561) q[0];
rz(-pi) q[1];
rz(1.4975411) q[2];
sx q[2];
rz(-3.0689291) q[2];
sx q[2];
rz(-2.6469028) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1637472) q[1];
sx q[1];
rz(-2.6085269) q[1];
sx q[1];
rz(-1.0372437) q[1];
x q[2];
rz(-2.1499885) q[3];
sx q[3];
rz(-2.0957855) q[3];
sx q[3];
rz(-0.92574173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8757214) q[2];
sx q[2];
rz(-3.044812) q[2];
sx q[2];
rz(-1.9270012) q[2];
rz(0.38861361) q[3];
sx q[3];
rz(-0.92258421) q[3];
sx q[3];
rz(-2.0428366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6471967) q[0];
sx q[0];
rz(-1.5460486) q[0];
sx q[0];
rz(-1.8478951) q[0];
rz(2.9248059) q[1];
sx q[1];
rz(-0.68991359) q[1];
sx q[1];
rz(0.51295113) q[1];
rz(-0.67736652) q[2];
sx q[2];
rz(-1.9147041) q[2];
sx q[2];
rz(-0.83821581) q[2];
rz(1.8372336) q[3];
sx q[3];
rz(-1.1188917) q[3];
sx q[3];
rz(1.6592142) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
