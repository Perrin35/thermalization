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
rz(2.9945381) q[0];
sx q[0];
rz(-1.7400063) q[0];
sx q[0];
rz(2.2214878) q[0];
rz(3.0609581) q[1];
sx q[1];
rz(-0.57890761) q[1];
sx q[1];
rz(-0.97715598) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29229627) q[0];
sx q[0];
rz(-1.9417282) q[0];
sx q[0];
rz(-1.3440668) q[0];
rz(3.009367) q[2];
sx q[2];
rz(-1.5572845) q[2];
sx q[2];
rz(0.15541645) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9639114) q[1];
sx q[1];
rz(-0.61574725) q[1];
sx q[1];
rz(-2.9065688) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8771873) q[3];
sx q[3];
rz(-0.9054817) q[3];
sx q[3];
rz(0.09395919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.5620293) q[2];
sx q[2];
rz(-1.2144438) q[2];
sx q[2];
rz(-0.62082949) q[2];
rz(0.51554716) q[3];
sx q[3];
rz(-1.4225057) q[3];
sx q[3];
rz(-0.8425042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8949378) q[0];
sx q[0];
rz(-1.6064914) q[0];
sx q[0];
rz(-0.57902336) q[0];
rz(0.82194263) q[1];
sx q[1];
rz(-1.6252981) q[1];
sx q[1];
rz(-2.6878405) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653719) q[0];
sx q[0];
rz(-2.1848618) q[0];
sx q[0];
rz(-1.8472415) q[0];
rz(-pi) q[1];
rz(0.2208925) q[2];
sx q[2];
rz(-2.6651987) q[2];
sx q[2];
rz(1.9770789) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0099185506) q[1];
sx q[1];
rz(-1.446053) q[1];
sx q[1];
rz(-2.591955) q[1];
rz(-pi) q[2];
rz(-2.6964397) q[3];
sx q[3];
rz(-1.0694519) q[3];
sx q[3];
rz(-1.8568764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4072676) q[2];
sx q[2];
rz(-3.0027323) q[2];
sx q[2];
rz(-0.56488758) q[2];
rz(-0.35483739) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(-1.423665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9722209) q[0];
sx q[0];
rz(-0.2114978) q[0];
sx q[0];
rz(-2.5900904) q[0];
rz(0.36477271) q[1];
sx q[1];
rz(-0.33939895) q[1];
sx q[1];
rz(-0.50484467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0612978) q[0];
sx q[0];
rz(-2.7027931) q[0];
sx q[0];
rz(-0.3831692) q[0];
rz(-0.41270035) q[2];
sx q[2];
rz(-1.4849714) q[2];
sx q[2];
rz(2.7049261) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6712196) q[1];
sx q[1];
rz(-0.80802901) q[1];
sx q[1];
rz(-0.68617448) q[1];
rz(1.9471719) q[3];
sx q[3];
rz(-2.0224114) q[3];
sx q[3];
rz(2.268689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43368936) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(0.047253963) q[2];
rz(-2.3047678) q[3];
sx q[3];
rz(-0.72332007) q[3];
sx q[3];
rz(1.0251454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4243917) q[0];
sx q[0];
rz(-1.0294788) q[0];
sx q[0];
rz(-1.3264054) q[0];
rz(-1.4855509) q[1];
sx q[1];
rz(-1.0842208) q[1];
sx q[1];
rz(0.63492376) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50872749) q[0];
sx q[0];
rz(-1.5888402) q[0];
sx q[0];
rz(-0.0012854171) q[0];
x q[1];
rz(1.0718173) q[2];
sx q[2];
rz(-1.8868251) q[2];
sx q[2];
rz(1.8763148) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.30764515) q[1];
sx q[1];
rz(-1.9059423) q[1];
sx q[1];
rz(-0.076471523) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4621099) q[3];
sx q[3];
rz(-0.75089291) q[3];
sx q[3];
rz(-0.39284947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2403468) q[2];
sx q[2];
rz(-2.2525658) q[2];
sx q[2];
rz(1.3725613) q[2];
rz(-1.2076123) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(2.8670368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.32886252) q[0];
sx q[0];
rz(-0.48506081) q[0];
sx q[0];
rz(0.10511705) q[0];
rz(0.33991995) q[1];
sx q[1];
rz(-2.2272019) q[1];
sx q[1];
rz(1.0629268) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2824616) q[0];
sx q[0];
rz(-2.4300008) q[0];
sx q[0];
rz(-3.1200073) q[0];
x q[1];
rz(-1.8196443) q[2];
sx q[2];
rz(-2.3169059) q[2];
sx q[2];
rz(0.4363554) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8430994) q[1];
sx q[1];
rz(-1.5122422) q[1];
sx q[1];
rz(-1.5850409) q[1];
rz(-pi) q[2];
rz(2.0696409) q[3];
sx q[3];
rz(-1.5728647) q[3];
sx q[3];
rz(-0.75132124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0792599) q[2];
sx q[2];
rz(-1.3881114) q[2];
sx q[2];
rz(-1.3366535) q[2];
rz(-0.054232728) q[3];
sx q[3];
rz(-1.1905866) q[3];
sx q[3];
rz(-0.59468734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.9688251) q[0];
sx q[0];
rz(-2.4496138) q[0];
sx q[0];
rz(-1.927595) q[0];
rz(-1.0460151) q[1];
sx q[1];
rz(-2.0966625) q[1];
sx q[1];
rz(-0.85375839) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3873685) q[0];
sx q[0];
rz(-1.4367391) q[0];
sx q[0];
rz(3.0813974) q[0];
x q[1];
rz(0.98967239) q[2];
sx q[2];
rz(-2.2377491) q[2];
sx q[2];
rz(-3.0158531) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.56928192) q[1];
sx q[1];
rz(-1.0195059) q[1];
sx q[1];
rz(2.2612919) q[1];
rz(-pi) q[2];
rz(2.7925909) q[3];
sx q[3];
rz(-1.4274538) q[3];
sx q[3];
rz(0.81763148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13009109) q[2];
sx q[2];
rz(-1.4902196) q[2];
sx q[2];
rz(0.74756527) q[2];
rz(2.2085564) q[3];
sx q[3];
rz(-1.3676164) q[3];
sx q[3];
rz(-0.30195495) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1099243) q[0];
sx q[0];
rz(-2.1878991) q[0];
sx q[0];
rz(-0.64144301) q[0];
rz(2.172442) q[1];
sx q[1];
rz(-1.0698003) q[1];
sx q[1];
rz(2.0279121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3308709) q[0];
sx q[0];
rz(-1.4735876) q[0];
sx q[0];
rz(1.9770798) q[0];
x q[1];
rz(-0.93514438) q[2];
sx q[2];
rz(-2.0144267) q[2];
sx q[2];
rz(-2.005524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.42150527) q[1];
sx q[1];
rz(-1.0522138) q[1];
sx q[1];
rz(0.20629136) q[1];
rz(2.6212531) q[3];
sx q[3];
rz(-1.5410149) q[3];
sx q[3];
rz(-0.34572085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5668737) q[2];
sx q[2];
rz(-2.4428664) q[2];
sx q[2];
rz(-2.3835772) q[2];
rz(0.30592439) q[3];
sx q[3];
rz(-2.7829792) q[3];
sx q[3];
rz(2.6942286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7253983) q[0];
sx q[0];
rz(-0.60878009) q[0];
sx q[0];
rz(2.388227) q[0];
rz(2.2196409) q[1];
sx q[1];
rz(-1.7148858) q[1];
sx q[1];
rz(3.0070378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050722402) q[0];
sx q[0];
rz(-2.7742552) q[0];
sx q[0];
rz(-1.3274756) q[0];
rz(-pi) q[1];
rz(-2.806611) q[2];
sx q[2];
rz(-0.6689531) q[2];
sx q[2];
rz(-1.9425336) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0363732) q[1];
sx q[1];
rz(-1.4810307) q[1];
sx q[1];
rz(-1.5898934) q[1];
rz(-pi) q[2];
rz(0.67340019) q[3];
sx q[3];
rz(-1.2030501) q[3];
sx q[3];
rz(2.0665702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.77384633) q[2];
sx q[2];
rz(-1.6951122) q[2];
sx q[2];
rz(1.9473677) q[2];
rz(1.9959244) q[3];
sx q[3];
rz(-2.001389) q[3];
sx q[3];
rz(1.0660508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0399748) q[0];
sx q[0];
rz(-0.45926738) q[0];
sx q[0];
rz(2.7591144) q[0];
rz(2.3756012) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(-1.8006178) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95582286) q[0];
sx q[0];
rz(-0.18322028) q[0];
sx q[0];
rz(-0.51143666) q[0];
x q[1];
rz(-3.0485504) q[2];
sx q[2];
rz(-0.75936717) q[2];
sx q[2];
rz(-1.134128) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0611524) q[1];
sx q[1];
rz(-0.78394475) q[1];
sx q[1];
rz(2.817201) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3519456) q[3];
sx q[3];
rz(-0.36389458) q[3];
sx q[3];
rz(-2.4879251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0535023) q[2];
sx q[2];
rz(-0.92038766) q[2];
sx q[2];
rz(-3.0255393) q[2];
rz(-1.2608438) q[3];
sx q[3];
rz(-2.0033658) q[3];
sx q[3];
rz(1.6837696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5905404) q[0];
sx q[0];
rz(-0.11369471) q[0];
sx q[0];
rz(0.1846479) q[0];
rz(-0.91839904) q[1];
sx q[1];
rz(-1.5469488) q[1];
sx q[1];
rz(-1.437423) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6708095) q[0];
sx q[0];
rz(-0.52007404) q[0];
sx q[0];
rz(1.1519679) q[0];
rz(1.499648) q[2];
sx q[2];
rz(-0.84991036) q[2];
sx q[2];
rz(-2.3325553) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5823707) q[1];
sx q[1];
rz(-2.5032024) q[1];
sx q[1];
rz(1.8965782) q[1];
rz(-pi) q[2];
rz(-1.2380794) q[3];
sx q[3];
rz(-2.2196688) q[3];
sx q[3];
rz(2.1391265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.664428) q[2];
sx q[2];
rz(-1.8958586) q[2];
sx q[2];
rz(-2.5028382) q[2];
rz(0.81513682) q[3];
sx q[3];
rz(-0.4709979) q[3];
sx q[3];
rz(-2.7208929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7246134) q[0];
sx q[0];
rz(-1.2495578) q[0];
sx q[0];
rz(0.15171262) q[0];
rz(-2.3607415) q[1];
sx q[1];
rz(-1.6897222) q[1];
sx q[1];
rz(2.2471468) q[1];
rz(1.8506321) q[2];
sx q[2];
rz(-1.6879514) q[2];
sx q[2];
rz(0.25456706) q[2];
rz(-2.9853447) q[3];
sx q[3];
rz(-2.8946946) q[3];
sx q[3];
rz(1.3794086) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
