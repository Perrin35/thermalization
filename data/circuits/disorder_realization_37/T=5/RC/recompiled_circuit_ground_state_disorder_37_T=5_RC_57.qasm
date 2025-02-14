OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6310298) q[0];
sx q[0];
rz(3.6462311) q[0];
sx q[0];
rz(12.991821) q[0];
rz(0.57234859) q[1];
sx q[1];
rz(4.2387716) q[1];
sx q[1];
rz(7.4833202) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10353032) q[0];
sx q[0];
rz(-1.5611042) q[0];
sx q[0];
rz(-2.8907625) q[0];
x q[1];
rz(-2.8699371) q[2];
sx q[2];
rz(-1.9052398) q[2];
sx q[2];
rz(-1.7890499) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.26781005) q[1];
sx q[1];
rz(-1.3569731) q[1];
sx q[1];
rz(-0.38114433) q[1];
rz(2.4704504) q[3];
sx q[3];
rz(-1.6832388) q[3];
sx q[3];
rz(2.6813316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1537062) q[2];
sx q[2];
rz(-1.989863) q[2];
sx q[2];
rz(-0.64725867) q[2];
rz(0.17793947) q[3];
sx q[3];
rz(-1.2578332) q[3];
sx q[3];
rz(1.2037163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4030289) q[0];
sx q[0];
rz(-3.0572427) q[0];
sx q[0];
rz(-2.6179598) q[0];
rz(-1.9117484) q[1];
sx q[1];
rz(-2.3975027) q[1];
sx q[1];
rz(-2.8935166) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46892525) q[0];
sx q[0];
rz(-1.6771731) q[0];
sx q[0];
rz(2.399978) q[0];
x q[1];
rz(1.1311943) q[2];
sx q[2];
rz(-1.9705551) q[2];
sx q[2];
rz(-0.93610379) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9776002) q[1];
sx q[1];
rz(-1.6095785) q[1];
sx q[1];
rz(-1.7165403) q[1];
x q[2];
rz(0.66623293) q[3];
sx q[3];
rz(-1.1403475) q[3];
sx q[3];
rz(-2.7048049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5292042) q[2];
sx q[2];
rz(-2.2233621) q[2];
sx q[2];
rz(2.6017792) q[2];
rz(-0.16048935) q[3];
sx q[3];
rz(-1.5925946) q[3];
sx q[3];
rz(0.81015712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4983343) q[0];
sx q[0];
rz(-2.46471) q[0];
sx q[0];
rz(-3.0078889) q[0];
rz(-1.5191822) q[1];
sx q[1];
rz(-1.2956053) q[1];
sx q[1];
rz(0.79230961) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80817079) q[0];
sx q[0];
rz(-1.9229888) q[0];
sx q[0];
rz(2.9461224) q[0];
x q[1];
rz(-2.9970005) q[2];
sx q[2];
rz(-1.4404313) q[2];
sx q[2];
rz(-3.0693288) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6349259) q[1];
sx q[1];
rz(-2.185195) q[1];
sx q[1];
rz(-0.8393112) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30012975) q[3];
sx q[3];
rz(-1.5709166) q[3];
sx q[3];
rz(-0.67226582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.85492674) q[2];
sx q[2];
rz(-2.4714405) q[2];
sx q[2];
rz(2.3112042) q[2];
rz(-1.3487799) q[3];
sx q[3];
rz(-1.034779) q[3];
sx q[3];
rz(-2.5428298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26375672) q[0];
sx q[0];
rz(-0.97665518) q[0];
sx q[0];
rz(0.42309716) q[0];
rz(-1.9844203) q[1];
sx q[1];
rz(-1.4856228) q[1];
sx q[1];
rz(-2.5194936) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34919993) q[0];
sx q[0];
rz(-2.3696941) q[0];
sx q[0];
rz(2.1436439) q[0];
x q[1];
rz(-2.9329002) q[2];
sx q[2];
rz(-1.1596173) q[2];
sx q[2];
rz(1.5087939) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2036977) q[1];
sx q[1];
rz(-2.7337498) q[1];
sx q[1];
rz(1.5289343) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.29633) q[3];
sx q[3];
rz(-1.5912531) q[3];
sx q[3];
rz(-1.4485698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28079924) q[2];
sx q[2];
rz(-1.7223027) q[2];
sx q[2];
rz(2.5234176) q[2];
rz(1.482796) q[3];
sx q[3];
rz(-1.7890472) q[3];
sx q[3];
rz(-1.5084722) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16267714) q[0];
sx q[0];
rz(-1.8122883) q[0];
sx q[0];
rz(-0.83258122) q[0];
rz(-2.816448) q[1];
sx q[1];
rz(-0.99601662) q[1];
sx q[1];
rz(-0.064373374) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58080855) q[0];
sx q[0];
rz(-0.34750313) q[0];
sx q[0];
rz(2.7827713) q[0];
rz(-pi) q[1];
rz(-2.5413248) q[2];
sx q[2];
rz(-1.1986599) q[2];
sx q[2];
rz(-1.063571) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0400914) q[1];
sx q[1];
rz(-2.1143997) q[1];
sx q[1];
rz(0.60740791) q[1];
x q[2];
rz(-0.35803087) q[3];
sx q[3];
rz(-0.97409407) q[3];
sx q[3];
rz(-2.1880414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23318204) q[2];
sx q[2];
rz(-2.5711214) q[2];
sx q[2];
rz(-1.849966) q[2];
rz(2.6455247) q[3];
sx q[3];
rz(-0.78947624) q[3];
sx q[3];
rz(1.9991649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9555776) q[0];
sx q[0];
rz(-0.29009524) q[0];
sx q[0];
rz(-0.41906038) q[0];
rz(0.31173197) q[1];
sx q[1];
rz(-1.5429976) q[1];
sx q[1];
rz(2.9980803) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3449609) q[0];
sx q[0];
rz(-1.2141742) q[0];
sx q[0];
rz(-2.659647) q[0];
rz(-pi) q[1];
rz(-0.70563282) q[2];
sx q[2];
rz(-2.7876283) q[2];
sx q[2];
rz(-1.7662545) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12998768) q[1];
sx q[1];
rz(-1.2530348) q[1];
sx q[1];
rz(0.77316534) q[1];
rz(-pi) q[2];
rz(-2.6426598) q[3];
sx q[3];
rz(-0.29424516) q[3];
sx q[3];
rz(0.1732451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.781337) q[2];
sx q[2];
rz(-2.2402253) q[2];
sx q[2];
rz(-1.1330913) q[2];
rz(-2.0356778) q[3];
sx q[3];
rz(-1.4191041) q[3];
sx q[3];
rz(0.47479409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.356242) q[0];
sx q[0];
rz(-0.79371912) q[0];
sx q[0];
rz(-1.6957977) q[0];
rz(-0.71714199) q[1];
sx q[1];
rz(-1.278806) q[1];
sx q[1];
rz(2.380127) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8728947) q[0];
sx q[0];
rz(-0.94572645) q[0];
sx q[0];
rz(1.3891267) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.285073) q[2];
sx q[2];
rz(-1.9437143) q[2];
sx q[2];
rz(-0.97717092) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95132212) q[1];
sx q[1];
rz(-0.99748625) q[1];
sx q[1];
rz(3.0828397) q[1];
rz(-pi) q[2];
rz(-1.7588571) q[3];
sx q[3];
rz(-1.7283963) q[3];
sx q[3];
rz(3.0732875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.36934272) q[2];
sx q[2];
rz(-2.6476634) q[2];
sx q[2];
rz(2.210145) q[2];
rz(-2.5233968) q[3];
sx q[3];
rz(-1.5074916) q[3];
sx q[3];
rz(-0.96562323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44929993) q[0];
sx q[0];
rz(-1.3626784) q[0];
sx q[0];
rz(-1.3344673) q[0];
rz(0.5131228) q[1];
sx q[1];
rz(-2.158439) q[1];
sx q[1];
rz(-1.9122874) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052027651) q[0];
sx q[0];
rz(-2.2917245) q[0];
sx q[0];
rz(2.5682279) q[0];
rz(-1.8304873) q[2];
sx q[2];
rz(-1.6605111) q[2];
sx q[2];
rz(-3.1083687) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.38459435) q[1];
sx q[1];
rz(-1.4149067) q[1];
sx q[1];
rz(-1.5499359) q[1];
rz(-pi) q[2];
rz(-1.5531814) q[3];
sx q[3];
rz(-0.9086787) q[3];
sx q[3];
rz(2.3603242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8339771) q[2];
sx q[2];
rz(-0.5270842) q[2];
sx q[2];
rz(-0.21044593) q[2];
rz(-3.0757507) q[3];
sx q[3];
rz(-2.6688771) q[3];
sx q[3];
rz(-2.3426447) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.598572) q[0];
sx q[0];
rz(-2.4871171) q[0];
sx q[0];
rz(1.2731113) q[0];
rz(-1.2741362) q[1];
sx q[1];
rz(-0.68398634) q[1];
sx q[1];
rz(-0.88409105) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1691723) q[0];
sx q[0];
rz(-0.61686777) q[0];
sx q[0];
rz(1.5703809) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7124518) q[2];
sx q[2];
rz(-0.56641266) q[2];
sx q[2];
rz(2.8265116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99779656) q[1];
sx q[1];
rz(-1.9854768) q[1];
sx q[1];
rz(2.8383377) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8330634) q[3];
sx q[3];
rz(-1.8179699) q[3];
sx q[3];
rz(0.17251523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.990443) q[2];
sx q[2];
rz(-0.98439011) q[2];
sx q[2];
rz(1.3746064) q[2];
rz(2.9511342) q[3];
sx q[3];
rz(-2.3951267) q[3];
sx q[3];
rz(-1.9577352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-0.16354887) q[0];
sx q[0];
rz(-1.6885641) q[0];
sx q[0];
rz(1.2646041) q[0];
rz(-3.0679852) q[1];
sx q[1];
rz(-1.0818447) q[1];
sx q[1];
rz(-0.61990613) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6524175) q[0];
sx q[0];
rz(-2.3803854) q[0];
sx q[0];
rz(1.8106145) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61087278) q[2];
sx q[2];
rz(-1.9985285) q[2];
sx q[2];
rz(1.7868702) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3504847) q[1];
sx q[1];
rz(-1.7813761) q[1];
sx q[1];
rz(-1.3665707) q[1];
rz(-pi) q[2];
rz(-2.6905493) q[3];
sx q[3];
rz(-2.6313734) q[3];
sx q[3];
rz(-3.1376145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8668883) q[2];
sx q[2];
rz(-2.4394749) q[2];
sx q[2];
rz(0.14599027) q[2];
rz(-0.22249666) q[3];
sx q[3];
rz(-0.22495088) q[3];
sx q[3];
rz(2.4116624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9324026) q[0];
sx q[0];
rz(-1.9575735) q[0];
sx q[0];
rz(0.14458543) q[0];
rz(1.9119541) q[1];
sx q[1];
rz(-1.3508136) q[1];
sx q[1];
rz(-1.2523686) q[1];
rz(-2.8242883) q[2];
sx q[2];
rz(-2.5226421) q[2];
sx q[2];
rz(0.68963827) q[2];
rz(3.1093507) q[3];
sx q[3];
rz(-1.3606291) q[3];
sx q[3];
rz(1.2442333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
