OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0260789) q[0];
sx q[0];
rz(-1.6576515) q[0];
sx q[0];
rz(-2.8154362) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(-1.7915373) q[1];
sx q[1];
rz(-1.5426558) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21539772) q[0];
sx q[0];
rz(-0.30456671) q[0];
sx q[0];
rz(-0.79415168) q[0];
x q[1];
rz(-1.8250263) q[2];
sx q[2];
rz(-2.1581274) q[2];
sx q[2];
rz(1.4884782) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1118288) q[1];
sx q[1];
rz(-1.3297237) q[1];
sx q[1];
rz(1.8087216) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5845675) q[3];
sx q[3];
rz(-0.78409401) q[3];
sx q[3];
rz(-0.23628326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2797543) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(2.2564783) q[2];
rz(2.4195813) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(-0.0074145934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9279813) q[0];
sx q[0];
rz(-0.95887029) q[0];
sx q[0];
rz(-2.0425178) q[0];
rz(2.4765769) q[1];
sx q[1];
rz(-1.4140833) q[1];
sx q[1];
rz(-2.2639993) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7859902) q[0];
sx q[0];
rz(-0.88760932) q[0];
sx q[0];
rz(-1.3366633) q[0];
x q[1];
rz(1.4363102) q[2];
sx q[2];
rz(-1.9324537) q[2];
sx q[2];
rz(-2.550617) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6192012) q[1];
sx q[1];
rz(-0.11208216) q[1];
sx q[1];
rz(-1.2680608) q[1];
rz(0.33438501) q[3];
sx q[3];
rz(-1.9352203) q[3];
sx q[3];
rz(0.78727608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7426804) q[2];
sx q[2];
rz(-1.30554) q[2];
sx q[2];
rz(-1.1304643) q[2];
rz(1.8418664) q[3];
sx q[3];
rz(-1.2239417) q[3];
sx q[3];
rz(-1.6931504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.995342) q[0];
sx q[0];
rz(-1.8595707) q[0];
sx q[0];
rz(-0.28999844) q[0];
rz(2.4747804) q[1];
sx q[1];
rz(-1.0338444) q[1];
sx q[1];
rz(-0.07382948) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2436115) q[0];
sx q[0];
rz(-1.5092761) q[0];
sx q[0];
rz(-3.0520526) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2697271) q[2];
sx q[2];
rz(-2.7060894) q[2];
sx q[2];
rz(1.0812024) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2245582) q[1];
sx q[1];
rz(-2.9552166) q[1];
sx q[1];
rz(1.264155) q[1];
rz(1.4521493) q[3];
sx q[3];
rz(-1.9199748) q[3];
sx q[3];
rz(-2.196764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6510216) q[2];
sx q[2];
rz(-1.9475503) q[2];
sx q[2];
rz(-2.4948965) q[2];
rz(2.0329287) q[3];
sx q[3];
rz(-0.78287786) q[3];
sx q[3];
rz(-1.1289319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17764238) q[0];
sx q[0];
rz(-0.17340604) q[0];
sx q[0];
rz(1.1886764) q[0];
rz(2.1229318) q[1];
sx q[1];
rz(-0.97266346) q[1];
sx q[1];
rz(-1.4368988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1543717) q[0];
sx q[0];
rz(-2.4043596) q[0];
sx q[0];
rz(3.0043976) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3827299) q[2];
sx q[2];
rz(-2.7042411) q[2];
sx q[2];
rz(0.26668374) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1632299) q[1];
sx q[1];
rz(-2.0320738) q[1];
sx q[1];
rz(-0.39050885) q[1];
rz(-0.8021637) q[3];
sx q[3];
rz(-2.9608316) q[3];
sx q[3];
rz(0.67644955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58549762) q[2];
sx q[2];
rz(-1.7079587) q[2];
sx q[2];
rz(1.1222703) q[2];
rz(-1.026011) q[3];
sx q[3];
rz(-0.75338537) q[3];
sx q[3];
rz(-0.99075738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68200237) q[0];
sx q[0];
rz(-2.2408709) q[0];
sx q[0];
rz(0.37297747) q[0];
rz(0.22398082) q[1];
sx q[1];
rz(-1.9517027) q[1];
sx q[1];
rz(-1.3164828) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017942863) q[0];
sx q[0];
rz(-1.1153478) q[0];
sx q[0];
rz(2.7148867) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8793648) q[2];
sx q[2];
rz(-1.3272459) q[2];
sx q[2];
rz(-1.3893931) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5455294) q[1];
sx q[1];
rz(-1.5585594) q[1];
sx q[1];
rz(1.5188602) q[1];
x q[2];
rz(-0.54721197) q[3];
sx q[3];
rz(-2.0378761) q[3];
sx q[3];
rz(-2.3576749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.10107723) q[2];
sx q[2];
rz(-2.310029) q[2];
sx q[2];
rz(-0.80580795) q[2];
rz(0.53330437) q[3];
sx q[3];
rz(-2.0093982) q[3];
sx q[3];
rz(-2.1300952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.39078113) q[0];
sx q[0];
rz(-1.3178786) q[0];
sx q[0];
rz(0.090963013) q[0];
rz(0.85995752) q[1];
sx q[1];
rz(-1.1227612) q[1];
sx q[1];
rz(1.8213173) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66619191) q[0];
sx q[0];
rz(-1.1328508) q[0];
sx q[0];
rz(1.9802718) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7485511) q[2];
sx q[2];
rz(-2.1338935) q[2];
sx q[2];
rz(0.69001889) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7652313) q[1];
sx q[1];
rz(-0.77502854) q[1];
sx q[1];
rz(-2.2750957) q[1];
rz(2.8517013) q[3];
sx q[3];
rz(-0.93512669) q[3];
sx q[3];
rz(1.0155201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0423353) q[2];
sx q[2];
rz(-0.96747413) q[2];
sx q[2];
rz(2.5406204) q[2];
rz(0.48505923) q[3];
sx q[3];
rz(-0.22189134) q[3];
sx q[3];
rz(1.4453567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8619974) q[0];
sx q[0];
rz(-1.1789362) q[0];
sx q[0];
rz(2.5860508) q[0];
rz(-3.1069966) q[1];
sx q[1];
rz(-0.75841537) q[1];
sx q[1];
rz(1.7506036) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44156528) q[0];
sx q[0];
rz(-1.0374984) q[0];
sx q[0];
rz(-1.9576859) q[0];
x q[1];
rz(2.7752635) q[2];
sx q[2];
rz(-1.3372256) q[2];
sx q[2];
rz(-2.4457268) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2512867) q[1];
sx q[1];
rz(-1.444) q[1];
sx q[1];
rz(-2.2433953) q[1];
rz(-pi) q[2];
rz(-1.4963385) q[3];
sx q[3];
rz(-2.3450608) q[3];
sx q[3];
rz(1.9394685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.81327072) q[2];
sx q[2];
rz(-0.44181028) q[2];
sx q[2];
rz(2.7461046) q[2];
rz(-1.288712) q[3];
sx q[3];
rz(-1.6059395) q[3];
sx q[3];
rz(-0.66974631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46273461) q[0];
sx q[0];
rz(-2.8023219) q[0];
sx q[0];
rz(1.6495552) q[0];
rz(-2.18816) q[1];
sx q[1];
rz(-2.0326734) q[1];
sx q[1];
rz(-1.7038201) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6839562) q[0];
sx q[0];
rz(-2.5914765) q[0];
sx q[0];
rz(1.9253299) q[0];
rz(-pi) q[1];
rz(-1.4346801) q[2];
sx q[2];
rz(-1.4021177) q[2];
sx q[2];
rz(-0.31751925) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21282141) q[1];
sx q[1];
rz(-1.0911687) q[1];
sx q[1];
rz(-0.70478435) q[1];
rz(-pi) q[2];
rz(-2.0182761) q[3];
sx q[3];
rz(-1.5527225) q[3];
sx q[3];
rz(-1.7759089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4961204) q[2];
sx q[2];
rz(-2.7068553) q[2];
sx q[2];
rz(0.17871857) q[2];
rz(-2.2802165) q[3];
sx q[3];
rz(-1.9390315) q[3];
sx q[3];
rz(-0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9361967) q[0];
sx q[0];
rz(-1.9367171) q[0];
sx q[0];
rz(-2.0478915) q[0];
rz(2.4049092) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(1.0587943) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98204389) q[0];
sx q[0];
rz(-1.0315572) q[0];
sx q[0];
rz(0.3190785) q[0];
rz(-pi) q[1];
rz(-2.7339897) q[2];
sx q[2];
rz(-1.1317481) q[2];
sx q[2];
rz(2.0087194) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4150548) q[1];
sx q[1];
rz(-1.6931567) q[1];
sx q[1];
rz(-1.9499669) q[1];
x q[2];
rz(2.6685647) q[3];
sx q[3];
rz(-1.2086476) q[3];
sx q[3];
rz(-2.2688697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.27292353) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(-1.6607364) q[2];
rz(2.7311834) q[3];
sx q[3];
rz(-1.7216262) q[3];
sx q[3];
rz(2.9836392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8695628) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(-2.8503382) q[0];
rz(-0.60925305) q[1];
sx q[1];
rz(-1.4657425) q[1];
sx q[1];
rz(1.7094918) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8490484) q[0];
sx q[0];
rz(-1.8869072) q[0];
sx q[0];
rz(-2.6228117) q[0];
rz(2.6860793) q[2];
sx q[2];
rz(-1.3580772) q[2];
sx q[2];
rz(-2.6754975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0309279) q[1];
sx q[1];
rz(-2.1222161) q[1];
sx q[1];
rz(1.5196147) q[1];
rz(-pi) q[2];
rz(-1.7581975) q[3];
sx q[3];
rz(-0.43863505) q[3];
sx q[3];
rz(2.8965829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6802406) q[2];
sx q[2];
rz(-1.1169008) q[2];
sx q[2];
rz(-2.0001901) q[2];
rz(1.6067778) q[3];
sx q[3];
rz(-1.9635868) q[3];
sx q[3];
rz(-0.46943584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5466945) q[0];
sx q[0];
rz(-1.5190769) q[0];
sx q[0];
rz(1.4357823) q[0];
rz(-2.7728511) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(-0.73451191) q[2];
sx q[2];
rz(-0.30167689) q[2];
sx q[2];
rz(-2.6405356) q[2];
rz(0.41714824) q[3];
sx q[3];
rz(-2.683831) q[3];
sx q[3];
rz(2.8575069) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
