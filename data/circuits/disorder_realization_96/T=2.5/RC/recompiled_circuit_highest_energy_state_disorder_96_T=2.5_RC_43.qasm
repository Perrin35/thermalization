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
rz(1.9198298) q[0];
sx q[0];
rz(-1.3718995) q[0];
sx q[0];
rz(-2.0576117) q[0];
rz(2.3594175) q[1];
sx q[1];
rz(-0.44936925) q[1];
sx q[1];
rz(-0.78701293) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1469022) q[0];
sx q[0];
rz(-2.7222607) q[0];
sx q[0];
rz(1.9492784) q[0];
x q[1];
rz(1.5151565) q[2];
sx q[2];
rz(-0.67650992) q[2];
sx q[2];
rz(1.8970053) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8690774) q[1];
sx q[1];
rz(-2.1109588) q[1];
sx q[1];
rz(0.16611986) q[1];
rz(-2.8008599) q[3];
sx q[3];
rz(-0.91487802) q[3];
sx q[3];
rz(2.4994196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7996173) q[2];
sx q[2];
rz(-1.3222597) q[2];
sx q[2];
rz(-2.4413617) q[2];
rz(-1.3618943) q[3];
sx q[3];
rz(-1.2075862) q[3];
sx q[3];
rz(0.12968682) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6561061) q[0];
sx q[0];
rz(-0.34657297) q[0];
sx q[0];
rz(1.19278) q[0];
rz(-2.9876409) q[1];
sx q[1];
rz(-0.62869453) q[1];
sx q[1];
rz(1.2492294) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.10478) q[0];
sx q[0];
rz(-1.6333155) q[0];
sx q[0];
rz(2.9680802) q[0];
x q[1];
rz(3.0871816) q[2];
sx q[2];
rz(-0.93428946) q[2];
sx q[2];
rz(-0.3776224) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8802048) q[1];
sx q[1];
rz(-1.6515762) q[1];
sx q[1];
rz(0.30552633) q[1];
x q[2];
rz(-0.36954255) q[3];
sx q[3];
rz(-2.2033354) q[3];
sx q[3];
rz(0.12128017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.93231702) q[2];
sx q[2];
rz(-1.2195769) q[2];
sx q[2];
rz(0.90927124) q[2];
rz(3.0026109) q[3];
sx q[3];
rz(-0.2551955) q[3];
sx q[3];
rz(2.2709258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22055498) q[0];
sx q[0];
rz(-1.166936) q[0];
sx q[0];
rz(-1.1880818) q[0];
rz(-2.2782169) q[1];
sx q[1];
rz(-1.8977576) q[1];
sx q[1];
rz(-0.9695425) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090102948) q[0];
sx q[0];
rz(-1.7307407) q[0];
sx q[0];
rz(-1.7311726) q[0];
rz(0.34805426) q[2];
sx q[2];
rz(-1.2173941) q[2];
sx q[2];
rz(-2.5329593) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.88501334) q[1];
sx q[1];
rz(-1.6672214) q[1];
sx q[1];
rz(0.23553358) q[1];
rz(2.2754955) q[3];
sx q[3];
rz(-1.1148165) q[3];
sx q[3];
rz(-0.34581071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7748176) q[2];
sx q[2];
rz(-2.884951) q[2];
sx q[2];
rz(0.65566629) q[2];
rz(1.5244779) q[3];
sx q[3];
rz(-1.3527801) q[3];
sx q[3];
rz(0.85339439) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26576385) q[0];
sx q[0];
rz(-1.1273071) q[0];
sx q[0];
rz(1.458459) q[0];
rz(-1.726285) q[1];
sx q[1];
rz(-1.4348607) q[1];
sx q[1];
rz(-1.8728135) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7527217) q[0];
sx q[0];
rz(-1.224504) q[0];
sx q[0];
rz(2.4779921) q[0];
x q[1];
rz(0.88246302) q[2];
sx q[2];
rz(-0.75661406) q[2];
sx q[2];
rz(-2.7992835) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6626265) q[1];
sx q[1];
rz(-0.72196805) q[1];
sx q[1];
rz(-1.8501227) q[1];
rz(-pi) q[2];
rz(-1.9400408) q[3];
sx q[3];
rz(-2.8974468) q[3];
sx q[3];
rz(-0.63318726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.65712523) q[2];
sx q[2];
rz(-1.3544269) q[2];
sx q[2];
rz(-1.8787059) q[2];
rz(0.2028939) q[3];
sx q[3];
rz(-2.1094567) q[3];
sx q[3];
rz(1.8111022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48388594) q[0];
sx q[0];
rz(-1.0840451) q[0];
sx q[0];
rz(-2.6935691) q[0];
rz(1.9810642) q[1];
sx q[1];
rz(-2.0612165) q[1];
sx q[1];
rz(2.0652658) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1370476) q[0];
sx q[0];
rz(-2.2820615) q[0];
sx q[0];
rz(0.28202951) q[0];
rz(-pi) q[1];
rz(2.9975909) q[2];
sx q[2];
rz(-2.8051734) q[2];
sx q[2];
rz(-0.1937723) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9256968) q[1];
sx q[1];
rz(-1.5091245) q[1];
sx q[1];
rz(-0.036176763) q[1];
rz(-0.95376261) q[3];
sx q[3];
rz(-1.5542398) q[3];
sx q[3];
rz(-1.7615033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6059604) q[2];
sx q[2];
rz(-1.863402) q[2];
sx q[2];
rz(0.34206259) q[2];
rz(-2.0096807) q[3];
sx q[3];
rz(-2.0703546) q[3];
sx q[3];
rz(-2.7560077) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1256506) q[0];
sx q[0];
rz(-0.60369879) q[0];
sx q[0];
rz(1.1018671) q[0];
rz(2.8562538) q[1];
sx q[1];
rz(-0.97831786) q[1];
sx q[1];
rz(0.91011059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0285595) q[0];
sx q[0];
rz(-2.9601149) q[0];
sx q[0];
rz(0.89361287) q[0];
rz(-pi) q[1];
rz(3.0481045) q[2];
sx q[2];
rz(-1.6851236) q[2];
sx q[2];
rz(1.0766407) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31426281) q[1];
sx q[1];
rz(-1.1276961) q[1];
sx q[1];
rz(-3.0783975) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15575445) q[3];
sx q[3];
rz(-1.7475583) q[3];
sx q[3];
rz(-0.246941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3204699) q[2];
sx q[2];
rz(-1.7197101) q[2];
sx q[2];
rz(1.0538496) q[2];
rz(-2.7779135) q[3];
sx q[3];
rz(-1.4854919) q[3];
sx q[3];
rz(0.53500879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66154552) q[0];
sx q[0];
rz(-2.1286025) q[0];
sx q[0];
rz(-2.9841828) q[0];
rz(2.5120381) q[1];
sx q[1];
rz(-1.5903571) q[1];
sx q[1];
rz(-3.0622283) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4422121) q[0];
sx q[0];
rz(-1.6621502) q[0];
sx q[0];
rz(-0.070789074) q[0];
rz(-pi) q[1];
rz(1.7479701) q[2];
sx q[2];
rz(-1.1821772) q[2];
sx q[2];
rz(1.3874229) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0473898) q[1];
sx q[1];
rz(-1.4416639) q[1];
sx q[1];
rz(0.49603059) q[1];
rz(-pi) q[2];
rz(-1.5299876) q[3];
sx q[3];
rz(-0.58678484) q[3];
sx q[3];
rz(-0.87848488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51084149) q[2];
sx q[2];
rz(-1.8264822) q[2];
sx q[2];
rz(-2.9691479) q[2];
rz(-0.6423966) q[3];
sx q[3];
rz(-2.1746217) q[3];
sx q[3];
rz(0.35014686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.073535) q[0];
sx q[0];
rz(-1.8386766) q[0];
sx q[0];
rz(-2.5750343) q[0];
rz(-2.9134275) q[1];
sx q[1];
rz(-1.6888432) q[1];
sx q[1];
rz(-1.8295005) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5540422) q[0];
sx q[0];
rz(-1.8186079) q[0];
sx q[0];
rz(-2.3378563) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.117339) q[2];
sx q[2];
rz(-0.84716958) q[2];
sx q[2];
rz(3.0229508) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.560558) q[1];
sx q[1];
rz(-1.0315064) q[1];
sx q[1];
rz(-0.59625397) q[1];
rz(-pi) q[2];
rz(-1.4890461) q[3];
sx q[3];
rz(-0.39021947) q[3];
sx q[3];
rz(1.8070328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.18099004) q[2];
sx q[2];
rz(-1.484551) q[2];
sx q[2];
rz(-2.9385938) q[2];
rz(0.8391884) q[3];
sx q[3];
rz(-0.33792096) q[3];
sx q[3];
rz(-0.85116974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39076552) q[0];
sx q[0];
rz(-0.44742328) q[0];
sx q[0];
rz(0.14697337) q[0];
rz(-0.43288639) q[1];
sx q[1];
rz(-1.9501481) q[1];
sx q[1];
rz(1.8870707) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1330652) q[0];
sx q[0];
rz(-2.0957895) q[0];
sx q[0];
rz(1.1675444) q[0];
rz(-pi) q[1];
rz(-1.3225088) q[2];
sx q[2];
rz(-1.5380368) q[2];
sx q[2];
rz(2.0078307) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.43707284) q[1];
sx q[1];
rz(-2.7531248) q[1];
sx q[1];
rz(2.0411496) q[1];
x q[2];
rz(0.12740429) q[3];
sx q[3];
rz(-0.16949305) q[3];
sx q[3];
rz(-2.3338846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7811084) q[2];
sx q[2];
rz(-2.1238748) q[2];
sx q[2];
rz(-0.54070365) q[2];
rz(0.30132076) q[3];
sx q[3];
rz(-0.40073985) q[3];
sx q[3];
rz(-1.2979243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28894579) q[0];
sx q[0];
rz(-1.9444554) q[0];
sx q[0];
rz(2.4135015) q[0];
rz(0.51746619) q[1];
sx q[1];
rz(-1.6518075) q[1];
sx q[1];
rz(0.99658406) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48595062) q[0];
sx q[0];
rz(-2.7962755) q[0];
sx q[0];
rz(1.4929857) q[0];
rz(-pi) q[1];
rz(2.9884981) q[2];
sx q[2];
rz(-2.5544807) q[2];
sx q[2];
rz(-0.73995164) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4415295) q[1];
sx q[1];
rz(-0.25694381) q[1];
sx q[1];
rz(1.5154626) q[1];
x q[2];
rz(0.18325652) q[3];
sx q[3];
rz(-1.6695287) q[3];
sx q[3];
rz(-1.2399286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1835798) q[2];
sx q[2];
rz(-1.5778342) q[2];
sx q[2];
rz(0.070240423) q[2];
rz(3.1349414) q[3];
sx q[3];
rz(-0.17156048) q[3];
sx q[3];
rz(-0.10449617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9329279) q[0];
sx q[0];
rz(-1.7345971) q[0];
sx q[0];
rz(-1.4165184) q[0];
rz(2.3729462) q[1];
sx q[1];
rz(-1.9192764) q[1];
sx q[1];
rz(2.8142014) q[1];
rz(-1.745426) q[2];
sx q[2];
rz(-3.027109) q[2];
sx q[2];
rz(-0.30828044) q[2];
rz(1.205659) q[3];
sx q[3];
rz(-1.8237999) q[3];
sx q[3];
rz(0.0013838125) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
