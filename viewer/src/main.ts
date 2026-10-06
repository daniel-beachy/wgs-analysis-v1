import { mount } from 'svelte';
import './app.css';
import App from './App.svelte';
import { startHeartbeat } from './lib/data';

mount(App, { target: document.getElementById('app')! });
startHeartbeat();
