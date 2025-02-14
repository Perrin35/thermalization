OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.12299744) q[0];
sx q[0];
rz(-1.7641492) q[0];
sx q[0];
rz(-2.7121845) q[0];
rz(1.6027066) q[1];
sx q[1];
rz(-1.4338926) q[1];
sx q[1];
rz(1.5418928) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7943952) q[0];
sx q[0];
rz(-1.7849677) q[0];
sx q[0];
rz(1.2783575) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5355655) q[2];
sx q[2];
rz(-1.4270947) q[2];
sx q[2];
rz(1.8409539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0059389) q[1];
sx q[1];
rz(-0.42891132) q[1];
sx q[1];
rz(0.072838293) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36088636) q[3];
sx q[3];
rz(-1.295422) q[3];
sx q[3];
rz(-2.9831246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8371007) q[2];
sx q[2];
rz(-1.9383177) q[2];
sx q[2];
rz(-2.4486747) q[2];
rz(2.9414226) q[3];
sx q[3];
rz(-0.19014159) q[3];
sx q[3];
rz(-2.0076803) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8189341) q[0];
sx q[0];
rz(-2.942473) q[0];
sx q[0];
rz(-2.0767427) q[0];
rz(-0.62243593) q[1];
sx q[1];
rz(-1.7935926) q[1];
sx q[1];
rz(-0.49076864) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3929185) q[0];
sx q[0];
rz(-0.025500925) q[0];
sx q[0];
rz(0.37125094) q[0];
rz(-1.347001) q[2];
sx q[2];
rz(-1.9836805) q[2];
sx q[2];
rz(-3.0996291) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.789114) q[1];
sx q[1];
rz(-1.6660457) q[1];
sx q[1];
rz(2.386341) q[1];
x q[2];
rz(-0.048100483) q[3];
sx q[3];
rz(-1.401618) q[3];
sx q[3];
rz(-0.98911232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6530767) q[2];
sx q[2];
rz(-1.642903) q[2];
sx q[2];
rz(-2.901279) q[2];
rz(2.6416685) q[3];
sx q[3];
rz(-2.5559055) q[3];
sx q[3];
rz(1.2691931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-1.2568473) q[0];
sx q[0];
rz(-2.186543) q[0];
sx q[0];
rz(0.98168674) q[0];
rz(2.6495972) q[1];
sx q[1];
rz(-2.056608) q[1];
sx q[1];
rz(1.5031776) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5971736) q[0];
sx q[0];
rz(-2.0453296) q[0];
sx q[0];
rz(0.76961343) q[0];
x q[1];
rz(3.0122224) q[2];
sx q[2];
rz(-0.40310848) q[2];
sx q[2];
rz(1.0327686) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.99276421) q[1];
sx q[1];
rz(-1.1849313) q[1];
sx q[1];
rz(2.7817315) q[1];
x q[2];
rz(1.4731367) q[3];
sx q[3];
rz(-1.5916232) q[3];
sx q[3];
rz(-1.5690924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.65695277) q[2];
sx q[2];
rz(-0.52565614) q[2];
sx q[2];
rz(-3.0204115) q[2];
rz(-2.1335404) q[3];
sx q[3];
rz(-2.0262148) q[3];
sx q[3];
rz(-2.0602267) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0693531) q[0];
sx q[0];
rz(-0.00088748137) q[0];
sx q[0];
rz(-2.6278611) q[0];
rz(-2.1836102) q[1];
sx q[1];
rz(-1.1424516) q[1];
sx q[1];
rz(-1.6987919) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49794562) q[0];
sx q[0];
rz(-1.4063324) q[0];
sx q[0];
rz(0.41467675) q[0];
rz(-pi) q[1];
rz(1.3495693) q[2];
sx q[2];
rz(-1.0509247) q[2];
sx q[2];
rz(-1.9398361) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.084667) q[1];
sx q[1];
rz(-2.5206893) q[1];
sx q[1];
rz(-1.7123182) q[1];
rz(0.6616627) q[3];
sx q[3];
rz(-1.6473624) q[3];
sx q[3];
rz(0.14684248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.55502597) q[2];
sx q[2];
rz(-0.96379605) q[2];
sx q[2];
rz(1.8761934) q[2];
rz(-1.966656) q[3];
sx q[3];
rz(-2.411071) q[3];
sx q[3];
rz(-1.1635715) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88605276) q[0];
sx q[0];
rz(-2.9386254) q[0];
sx q[0];
rz(1.3245921) q[0];
rz(-2.1728204) q[1];
sx q[1];
rz(-1.6561534) q[1];
sx q[1];
rz(-2.103215) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0540349) q[0];
sx q[0];
rz(-1.731858) q[0];
sx q[0];
rz(2.7557082) q[0];
rz(-pi) q[1];
rz(2.4244196) q[2];
sx q[2];
rz(-2.178827) q[2];
sx q[2];
rz(-3.0423683) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5165625) q[1];
sx q[1];
rz(-1.5296401) q[1];
sx q[1];
rz(2.7709318) q[1];
x q[2];
rz(2.0821003) q[3];
sx q[3];
rz(-0.97701525) q[3];
sx q[3];
rz(2.4932414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5021299) q[2];
sx q[2];
rz(-0.30439964) q[2];
sx q[2];
rz(0.89140618) q[2];
rz(-1.8799479) q[3];
sx q[3];
rz(-1.6311878) q[3];
sx q[3];
rz(0.6663028) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0797743) q[0];
sx q[0];
rz(-0.49538716) q[0];
sx q[0];
rz(-0.99639446) q[0];
rz(0.90006104) q[1];
sx q[1];
rz(-0.78949094) q[1];
sx q[1];
rz(0.17734227) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19781659) q[0];
sx q[0];
rz(-1.7223486) q[0];
sx q[0];
rz(3.1413548) q[0];
rz(-1.5992237) q[2];
sx q[2];
rz(-0.6387944) q[2];
sx q[2];
rz(2.6306689) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9070322) q[1];
sx q[1];
rz(-1.8563358) q[1];
sx q[1];
rz(-2.0689059) q[1];
rz(3.0445339) q[3];
sx q[3];
rz(-1.0639816) q[3];
sx q[3];
rz(2.7926796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8511054) q[2];
sx q[2];
rz(-1.9600211) q[2];
sx q[2];
rz(2.8372852) q[2];
rz(0.32637063) q[3];
sx q[3];
rz(-2.5984952) q[3];
sx q[3];
rz(-0.91199005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0684763) q[0];
sx q[0];
rz(-0.40962064) q[0];
sx q[0];
rz(2.2090744) q[0];
rz(0.069123507) q[1];
sx q[1];
rz(-0.2176452) q[1];
sx q[1];
rz(-2.1014012) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5297663) q[0];
sx q[0];
rz(-1.5965723) q[0];
sx q[0];
rz(2.3759205) q[0];
rz(-1.3142775) q[2];
sx q[2];
rz(-1.4019012) q[2];
sx q[2];
rz(1.0948563) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8833163) q[1];
sx q[1];
rz(-2.3280716) q[1];
sx q[1];
rz(1.3431647) q[1];
rz(-pi) q[2];
rz(-2.4916441) q[3];
sx q[3];
rz(-2.3293265) q[3];
sx q[3];
rz(-2.5358729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0065464) q[2];
sx q[2];
rz(-0.73837215) q[2];
sx q[2];
rz(-2.4086003) q[2];
rz(2.0586355) q[3];
sx q[3];
rz(-2.0854009) q[3];
sx q[3];
rz(2.1347031) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55423823) q[0];
sx q[0];
rz(-0.08389689) q[0];
sx q[0];
rz(2.6485637) q[0];
rz(-2.9250277) q[1];
sx q[1];
rz(-2.000587) q[1];
sx q[1];
rz(-2.1176178) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1676583) q[0];
sx q[0];
rz(-1.5827145) q[0];
sx q[0];
rz(0.59597909) q[0];
rz(-0.33228428) q[2];
sx q[2];
rz(-1.4632676) q[2];
sx q[2];
rz(-0.33974732) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1223381) q[1];
sx q[1];
rz(-2.755909) q[1];
sx q[1];
rz(0.38094873) q[1];
x q[2];
rz(0.92864236) q[3];
sx q[3];
rz(-1.3285884) q[3];
sx q[3];
rz(2.1845269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0563125) q[2];
sx q[2];
rz(-1.5093426) q[2];
sx q[2];
rz(2.4658266) q[2];
rz(-0.41302776) q[3];
sx q[3];
rz(-1.690381) q[3];
sx q[3];
rz(-2.5954424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5125047) q[0];
sx q[0];
rz(-2.4830723) q[0];
sx q[0];
rz(-0.034544695) q[0];
rz(-3.0870364) q[1];
sx q[1];
rz(-1.9787534) q[1];
sx q[1];
rz(-0.82130718) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4865882) q[0];
sx q[0];
rz(-1.4711709) q[0];
sx q[0];
rz(-0.57488802) q[0];
x q[1];
rz(3.05282) q[2];
sx q[2];
rz(-1.5615511) q[2];
sx q[2];
rz(-2.451596) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7682338) q[1];
sx q[1];
rz(-1.4211417) q[1];
sx q[1];
rz(-1.6878769) q[1];
x q[2];
rz(2.6590632) q[3];
sx q[3];
rz(-2.165613) q[3];
sx q[3];
rz(-1.9687259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7811232) q[2];
sx q[2];
rz(-1.0162153) q[2];
sx q[2];
rz(3.0086369) q[2];
rz(2.7305056) q[3];
sx q[3];
rz(-1.6229595) q[3];
sx q[3];
rz(2.0558004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0947615) q[0];
sx q[0];
rz(-0.60698858) q[0];
sx q[0];
rz(1.2588311) q[0];
rz(1.5065441) q[1];
sx q[1];
rz(-0.656744) q[1];
sx q[1];
rz(-0.81319317) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5933605) q[0];
sx q[0];
rz(-1.5357067) q[0];
sx q[0];
rz(-0.078087383) q[0];
rz(-pi) q[1];
rz(-1.7178463) q[2];
sx q[2];
rz(-0.67710256) q[2];
sx q[2];
rz(-0.33237095) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.53467864) q[1];
sx q[1];
rz(-0.77086222) q[1];
sx q[1];
rz(-2.4669564) q[1];
x q[2];
rz(1.7045295) q[3];
sx q[3];
rz(-1.1220891) q[3];
sx q[3];
rz(2.1577842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4930111) q[2];
sx q[2];
rz(-1.1375256) q[2];
sx q[2];
rz(1.0515593) q[2];
rz(-1.018853) q[3];
sx q[3];
rz(-2.2994883) q[3];
sx q[3];
rz(0.65783182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9857585) q[0];
sx q[0];
rz(-0.65407615) q[0];
sx q[0];
rz(0.88055897) q[0];
rz(-1.3195994) q[1];
sx q[1];
rz(-1.9965912) q[1];
sx q[1];
rz(0.051864787) q[1];
rz(-2.8116799) q[2];
sx q[2];
rz(-0.22116667) q[2];
sx q[2];
rz(0.55851182) q[2];
rz(0.29215688) q[3];
sx q[3];
rz(-2.9360129) q[3];
sx q[3];
rz(2.871411) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
