OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55460632) q[0];
sx q[0];
rz(2.245683) q[0];
sx q[0];
rz(10.829344) q[0];
rz(-2.4401234) q[1];
sx q[1];
rz(-2.520732) q[1];
sx q[1];
rz(-1.4863185) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1370927) q[0];
sx q[0];
rz(-0.21919964) q[0];
sx q[0];
rz(-0.57114925) q[0];
x q[1];
rz(1.4990184) q[2];
sx q[2];
rz(-1.992618) q[2];
sx q[2];
rz(-0.81411241) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9857786) q[1];
sx q[1];
rz(-1.474838) q[1];
sx q[1];
rz(2.7491261) q[1];
rz(0.078943723) q[3];
sx q[3];
rz(-1.5642089) q[3];
sx q[3];
rz(1.6062669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0861337) q[2];
sx q[2];
rz(-1.7666585) q[2];
sx q[2];
rz(1.3936183) q[2];
rz(2.3404775) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(1.0629517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.79008094) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(3.063607) q[0];
rz(2.4474735) q[1];
sx q[1];
rz(-2.0072939) q[1];
sx q[1];
rz(2.7819113) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168424) q[0];
sx q[0];
rz(-2.1152088) q[0];
sx q[0];
rz(2.398688) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7799136) q[2];
sx q[2];
rz(-2.2799727) q[2];
sx q[2];
rz(1.4271971) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2962131) q[1];
sx q[1];
rz(-1.737533) q[1];
sx q[1];
rz(2.1026033) q[1];
rz(-pi) q[2];
rz(-0.85675591) q[3];
sx q[3];
rz(-1.6008915) q[3];
sx q[3];
rz(-1.4853958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.89547196) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(0.94334156) q[2];
rz(1.881276) q[3];
sx q[3];
rz(-2.3337119) q[3];
sx q[3];
rz(-2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27186069) q[0];
sx q[0];
rz(-1.6534709) q[0];
sx q[0];
rz(0.46762064) q[0];
rz(0.46472654) q[1];
sx q[1];
rz(-0.48370353) q[1];
sx q[1];
rz(0.31276774) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62374672) q[0];
sx q[0];
rz(-0.22652921) q[0];
sx q[0];
rz(1.2073327) q[0];
rz(-1.9415641) q[2];
sx q[2];
rz(-2.3808378) q[2];
sx q[2];
rz(2.8460381) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86381972) q[1];
sx q[1];
rz(-1.7120655) q[1];
sx q[1];
rz(-1.3978676) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7346482) q[3];
sx q[3];
rz(-2.787628) q[3];
sx q[3];
rz(1.0696166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9096845) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(-2.8364733) q[2];
rz(-1.8049847) q[3];
sx q[3];
rz(-0.87953416) q[3];
sx q[3];
rz(0.46521503) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6453648) q[0];
sx q[0];
rz(-0.35571686) q[0];
sx q[0];
rz(-0.028045068) q[0];
rz(-1.4656981) q[1];
sx q[1];
rz(-0.60246712) q[1];
sx q[1];
rz(0.67273295) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4417277) q[0];
sx q[0];
rz(-1.8687975) q[0];
sx q[0];
rz(2.5900048) q[0];
rz(1.8276617) q[2];
sx q[2];
rz(-1.0360403) q[2];
sx q[2];
rz(1.9321439) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.67191254) q[1];
sx q[1];
rz(-1.8021936) q[1];
sx q[1];
rz(-3.0242821) q[1];
rz(-pi) q[2];
rz(-2.1898515) q[3];
sx q[3];
rz(-1.674106) q[3];
sx q[3];
rz(2.2593474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4611886) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(2.380774) q[2];
rz(-0.86756724) q[3];
sx q[3];
rz(-1.8224199) q[3];
sx q[3];
rz(0.76604617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5701533) q[0];
sx q[0];
rz(-2.0338991) q[0];
sx q[0];
rz(2.0779579) q[0];
rz(1.6617552) q[1];
sx q[1];
rz(-0.96264833) q[1];
sx q[1];
rz(-0.038287727) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83513658) q[0];
sx q[0];
rz(-1.4493754) q[0];
sx q[0];
rz(1.2760713) q[0];
x q[1];
rz(0.64735909) q[2];
sx q[2];
rz(-2.416496) q[2];
sx q[2];
rz(0.65078306) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8582279) q[1];
sx q[1];
rz(-1.3971395) q[1];
sx q[1];
rz(-3.0221992) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0376301) q[3];
sx q[3];
rz(-1.8133014) q[3];
sx q[3];
rz(-1.2638826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9866207) q[2];
sx q[2];
rz(-1.203323) q[2];
sx q[2];
rz(-1.0920452) q[2];
rz(1.6070222) q[3];
sx q[3];
rz(-0.97073308) q[3];
sx q[3];
rz(0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9606278) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(1.3076179) q[0];
rz(3.0788105) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(0.19097701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073478621) q[0];
sx q[0];
rz(-1.3604135) q[0];
sx q[0];
rz(0.14260261) q[0];
rz(-0.88626185) q[2];
sx q[2];
rz(-1.2949416) q[2];
sx q[2];
rz(-2.4133298) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.2497017) q[1];
sx q[1];
rz(-0.54348031) q[1];
sx q[1];
rz(1.5312503) q[1];
rz(-2.7950068) q[3];
sx q[3];
rz(-2.1086018) q[3];
sx q[3];
rz(-1.1188521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.442231) q[2];
sx q[2];
rz(-1.3012393) q[2];
sx q[2];
rz(2.7065275) q[2];
rz(-0.89138952) q[3];
sx q[3];
rz(-1.9458709) q[3];
sx q[3];
rz(1.58266) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94423914) q[0];
sx q[0];
rz(-0.22560142) q[0];
sx q[0];
rz(-0.42022589) q[0];
rz(-2.6603783) q[1];
sx q[1];
rz(-0.88164202) q[1];
sx q[1];
rz(1.0302672) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97354613) q[0];
sx q[0];
rz(-1.1991812) q[0];
sx q[0];
rz(-2.6269873) q[0];
x q[1];
rz(-1.9040362) q[2];
sx q[2];
rz(-0.2313279) q[2];
sx q[2];
rz(1.7817792) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8458) q[1];
sx q[1];
rz(-2.7880221) q[1];
sx q[1];
rz(-2.5938864) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3805192) q[3];
sx q[3];
rz(-1.7576302) q[3];
sx q[3];
rz(0.877617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.39945012) q[2];
sx q[2];
rz(-0.43562499) q[2];
sx q[2];
rz(-0.31526652) q[2];
rz(2.1049843) q[3];
sx q[3];
rz(-1.2549812) q[3];
sx q[3];
rz(2.172327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2675562) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(-1.9497005) q[0];
rz(-2.2299956) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(-2.079516) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5870283) q[0];
sx q[0];
rz(-1.4438859) q[0];
sx q[0];
rz(2.1112604) q[0];
rz(2.383963) q[2];
sx q[2];
rz(-0.55765753) q[2];
sx q[2];
rz(0.66495313) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8909104) q[1];
sx q[1];
rz(-0.85103304) q[1];
sx q[1];
rz(2.2834884) q[1];
x q[2];
rz(1.6154556) q[3];
sx q[3];
rz(-2.4756845) q[3];
sx q[3];
rz(2.8093616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9553817) q[2];
sx q[2];
rz(-0.28327981) q[2];
sx q[2];
rz(-0.051076802) q[2];
rz(-0.71930277) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1445769) q[0];
sx q[0];
rz(-2.7269195) q[0];
sx q[0];
rz(-0.64494079) q[0];
rz(-2.0058517) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(2.2924246) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5727974) q[0];
sx q[0];
rz(-2.704247) q[0];
sx q[0];
rz(2.1225147) q[0];
rz(-0.37306771) q[2];
sx q[2];
rz(-0.39857769) q[2];
sx q[2];
rz(-1.0363491) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2850212) q[1];
sx q[1];
rz(-1.3339431) q[1];
sx q[1];
rz(2.9076438) q[1];
x q[2];
rz(1.4056021) q[3];
sx q[3];
rz(-1.7690036) q[3];
sx q[3];
rz(2.7310731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7944305) q[2];
sx q[2];
rz(-0.55134761) q[2];
sx q[2];
rz(0.17865044) q[2];
rz(0.84154877) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(2.6508022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0586044) q[0];
sx q[0];
rz(-1.2036136) q[0];
sx q[0];
rz(2.2264746) q[0];
rz(1.2471584) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(2.9291005) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8520078) q[0];
sx q[0];
rz(-0.6913018) q[0];
sx q[0];
rz(1.9612802) q[0];
x q[1];
rz(3.100527) q[2];
sx q[2];
rz(-1.8827056) q[2];
sx q[2];
rz(1.8598156) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5455268) q[1];
sx q[1];
rz(-2.3061228) q[1];
sx q[1];
rz(2.2722785) q[1];
rz(-pi) q[2];
rz(-0.28189567) q[3];
sx q[3];
rz(-2.3882139) q[3];
sx q[3];
rz(0.58338141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.96182573) q[2];
sx q[2];
rz(-0.13040725) q[2];
sx q[2];
rz(-1.4546222) q[2];
rz(-1.9466594) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(2.6132244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4029978) q[0];
sx q[0];
rz(-1.4017372) q[0];
sx q[0];
rz(-1.5177939) q[0];
rz(0.075642792) q[1];
sx q[1];
rz(-1.5374001) q[1];
sx q[1];
rz(-1.7061445) q[1];
rz(-3.1235789) q[2];
sx q[2];
rz(-1.6406888) q[2];
sx q[2];
rz(-2.3802118) q[2];
rz(-2.6247737) q[3];
sx q[3];
rz(-1.9782981) q[3];
sx q[3];
rz(0.78936418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
